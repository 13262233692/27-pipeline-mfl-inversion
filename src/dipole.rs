use rayon::prelude::*;
use crossbeam_channel::{Sender, bounded};
use std::sync::atomic::{AtomicU64, Ordering};

use crate::types::{
    DefectProfile, DefectSeverity, PipeScanSegment, WALL_THICKNESS_M,
    MAGNETIC_PERMEABILITY,
};

const LIFT_OFF_M: f64 = 0.025;
const DIPOLE_GRID_SIZE: usize = 21;
const INVERSION_MAX_ITER: usize = 100;
const INVERSION_TOLERANCE: f64 = 1e-8;
const GAUSSIAN_FILTER_WINDOW: usize = 5;
const MAGNETIZATION_A_M: f64 = 1.0e6;
const NOISE_THRESHOLD_T: f64 = 5.0e-6;

const PROGRESS_THROTTLE_EVERY: u64 = 256;

pub enum ProgressMsg {
    StageReading(f64),
    StageInverting(f64),
    StageStats(f64),
    Done,
}

pub struct ProgressReporter {
    sender: Sender<ProgressMsg>,
    invert_counter: AtomicU64,
    read_counter: AtomicU64,
    total_work: u64,
}

impl ProgressReporter {
    pub fn new(sender: Sender<ProgressMsg>, total_work: u64) -> Self {
        Self {
            sender,
            invert_counter: AtomicU64::new(0),
            read_counter: AtomicU64::new(0),
            total_work,
        }
    }

    #[inline]
    pub fn tick_read(&self, delta: u64) {
        let prev = self.read_counter.fetch_add(delta, Ordering::Relaxed);
        if (prev + delta) / PROGRESS_THROTTLE_EVERY != prev / PROGRESS_THROTTLE_EVERY {
            let ratio = (prev + delta) as f64 / self.total_work.max(1) as f64;
            let _ = self.sender.try_send(ProgressMsg::StageReading(ratio.min(1.0)));
        }
    }

    #[inline]
    pub fn tick_invert(&self, delta: u64) {
        let prev = self.invert_counter.fetch_add(delta, Ordering::Relaxed);
        if (prev + delta) / PROGRESS_THROTTLE_EVERY != prev / PROGRESS_THROTTLE_EVERY {
            let ratio = (prev + delta) as f64 / self.total_work.max(1) as f64;
            let _ = self.sender.try_send(ProgressMsg::StageInverting(ratio.min(1.0)));
        }
    }

    pub fn done(&self) {
        let _ = self.sender.send(ProgressMsg::Done);
    }
}

pub fn create_progress_channel() -> (Sender<ProgressMsg>, crossbeam_channel::Receiver<ProgressMsg>) {
    bounded::<ProgressMsg>(16)
}

pub struct DipoleInverter {
    wall_thickness: f64,
    lift_off: f64,
    permeability: f64,
    max_iterations: usize,
    tolerance: f64,
}

impl Default for DipoleInverter {
    fn default() -> Self {
        Self {
            wall_thickness: WALL_THICKNESS_M,
            lift_off: LIFT_OFF_M,
            permeability: MAGNETIC_PERMEABILITY,
            max_iterations: INVERSION_MAX_ITER,
            tolerance: INVERSION_TOLERANCE,
        }
    }
}

impl DipoleInverter {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_wall_thickness(mut self, thickness: f64) -> Self {
        self.wall_thickness = thickness;
        self
    }

    pub fn invert_segment_parallel(&self, segment: &mut PipeScanSegment) {
        self.invert_segment_parallel_with_progress(segment, None)
    }

    pub fn invert_segment_parallel_with_progress(
        &self,
        segment: &mut PipeScanSegment,
        reporter: Option<&ProgressReporter>,
    ) {
        let num_sensors = segment.sensor_data.len();
        let num_samples = segment.sensor_data[0].len();

        let mut defect_map: Vec<Vec<DefectProfile>> = vec![vec![DefectProfile::default(); num_samples]; num_sensors];

        defect_map.par_iter_mut().enumerate().for_each(|(sensor_idx, sensor_row)| {
            let signal: Vec<f64> = (0..num_samples)
                .map(|s| segment.sensor_data[sensor_idx][s].radial)
                .collect();

            let baseline = estimate_baseline(&signal);
            let anomaly: Vec<f64> = signal.iter().map(|&v| v - baseline).collect();
            let filtered = gaussian_filter(&anomaly, GAUSSIAN_FILTER_WINDOW);

            let mut local_ticks: u64 = 0;

            for s in 0..num_samples {
                let anomaly_magnitude = filtered[s].abs();

                let (depth, axial_length) = self.invert_single_point(
                    anomaly_magnitude,
                    &filtered,
                    s,
                    sensor_idx,
                    segment.axial_resolution_m(),
                );

                let severity = DefectSeverity::from_depth_ratio(depth / self.wall_thickness);

                sensor_row[s] = DefectProfile {
                    depth_m: depth,
                    axial_length_m: axial_length,
                    severity,
                };

                local_ticks += 1;
                if local_ticks >= PROGRESS_THROTTLE_EVERY {
                    if let Some(r) = reporter {
                        r.tick_invert(local_ticks);
                    }
                    local_ticks = 0;
                }
            }

            if local_ticks > 0 {
                if let Some(r) = reporter {
                    r.tick_invert(local_ticks);
                }
            }
        });

        segment.defect_map = defect_map;
    }

    fn invert_single_point(
        &self,
        anomaly_magnitude: f64,
        signal: &[f64],
        center_idx: usize,
        _sensor_idx: usize,
        axial_resolution: f64,
    ) -> (f64, f64) {
        if anomaly_magnitude < NOISE_THRESHOLD_T || !anomaly_magnitude.is_finite() {
            return (0.0, 0.0);
        }

        let safe_center = center_idx.min(signal.len().saturating_sub(1));

        let axial_length = estimate_axial_length(signal, safe_center, axial_resolution);

        if !axial_length.is_finite() || axial_length <= 0.0 {
            return (0.0, 0.0);
        }

        let depth = self.dipole_depth_inversion(anomaly_magnitude, axial_length);

        if depth.is_finite() && depth > 0.0 {
            (depth.min(self.wall_thickness * 0.95).max(0.0), axial_length)
        } else {
            (0.0, 0.0)
        }
    }

    fn dipole_depth_inversion(&self, b_radial_peak: f64, axial_length: f64) -> f64 {
        if b_radial_peak < 1e-12 || axial_length < 1e-6 {
            return 0.0;
        }
        if !b_radial_peak.is_finite() || !axial_length.is_finite() {
            return 0.0;
        }

        let target = b_radial_peak.abs();

        let b_max = self.forward_dipole_field(self.wall_thickness * 0.99, axial_length);
        if !b_max.is_finite() || target >= b_max {
            return self.wall_thickness * 0.99;
        }

        let b_min = self.forward_dipole_field(self.wall_thickness * 0.01, axial_length);
        if !b_min.is_finite() || target <= b_min {
            return self.wall_thickness * 0.01;
        }

        if b_max <= b_min {
            return self.wall_thickness * 0.5;
        }

        let mut low = self.wall_thickness * 0.01;
        let mut high = self.wall_thickness * 0.99;

        let mut iter_count = 0;
        while iter_count < self.max_iterations {
            let mid = (low + high) * 0.5;
            let b_mid = self.forward_dipole_field(mid, axial_length);

            if !b_mid.is_finite() {
                break;
            }

            let rel_err = (b_mid - target).abs() / target;
            if rel_err < self.tolerance {
                return mid;
            }

            if b_mid < target {
                low = mid;
            } else {
                high = mid;
            }

            if (high - low) / self.wall_thickness < 1e-6 {
                break;
            }
            iter_count += 1;
        }

        (low + high) * 0.5
    }

    fn forward_dipole_field(&self, depth: f64, axial_length: f64) -> f64 {
        if !depth.is_finite() || !axial_length.is_finite() {
            return 0.0;
        }

        let safe_depth = depth.max(1e-9).min(self.wall_thickness);
        let safe_axial = axial_length.max(1e-9);

        let z_sensor = self.wall_thickness + self.lift_off;
        let z_defect_top = self.wall_thickness - safe_depth;
        let z_center = z_defect_top + safe_depth / 2.0;

        let dz = z_sensor - z_center;
        let half_length = safe_axial / 2.0;

        let num_dipoles = DIPOLE_GRID_SIZE;
        let mut total_field = 0.0;
        let inv_4pi = self.permeability / (4.0 * std::f64::consts::PI);
        let m_z = safe_depth * 2.0;

        for i in 0..num_dipoles {
            let frac = (i as f64 + 0.5) / num_dipoles as f64;
            let x_dipole = -half_length + frac * safe_axial;

            let r_sq = x_dipole * x_dipole + dz * dz;
            if r_sq < 1e-18 {
                continue;
            }
            let r = r_sq.sqrt();
            let r_5 = r_sq * r_sq * r;
            if r_5 < 1e-45 {
                continue;
            }

            let term_3dz2_rsq = 3.0 * dz * dz - r_sq;

            let b_radial = inv_4pi * m_z * term_3dz2_rsq / r_5;

            if b_radial.is_finite() {
                total_field += b_radial;
            }
        }

        let avg = total_field / num_dipoles as f64;
        if avg.is_finite() { avg } else { 0.0 }
    }
}

fn estimate_baseline(signal: &[f64]) -> f64 {
    if signal.is_empty() {
        return 0.0;
    }

    let mut cleaned: Vec<f64> = signal
        .iter()
        .copied()
        .filter(|&v| v.is_finite())
        .collect();

    if cleaned.is_empty() {
        return 0.0;
    }

    cleaned.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let n = cleaned.len();
    let lower_quartile = cleaned[n / 4];
    let upper_quartile = cleaned[3 * n / 4];
    let iqr = upper_quartile - lower_quartile;

    if !iqr.is_finite() || iqr.abs() < 1e-18 {
        return cleaned[n / 2];
    }

    let filtered: Vec<f64> = cleaned
        .iter()
        .filter(|&&v| v >= lower_quartile - 1.5 * iqr && v <= upper_quartile + 1.5 * iqr)
        .copied()
        .collect();

    if filtered.is_empty() {
        return cleaned[n / 2];
    }

    filtered.iter().sum::<f64>() / filtered.len() as f64
}

fn gaussian_filter(signal: &[f64], window_size: usize) -> Vec<f64> {
    let n = signal.len();
    if n == 0 {
        return Vec::new();
    }

    let half = window_size / 2;
    let mut kernel = vec![0.0f64; window_size];
    let sigma = window_size as f64 / 6.0;

    for i in 0..window_size {
        let x = (i as f64 - half as f64) / sigma;
        kernel[i] = (-0.5 * x * x).exp();
    }

    let sum: f64 = kernel.iter().sum();
    if sum.abs() < 1e-18 {
        return signal.to_vec();
    }
    for k in &mut kernel {
        *k /= sum;
    }

    let mut result = vec![0.0f64; n];

    for i in 0..n {
        let mut val = 0.0;
        for j in 0..window_size {
            let idx = (i as isize - half as isize + j as isize).max(0).min(n as isize - 1) as usize;
            let sv = signal[idx];
            if sv.is_finite() {
                val += sv * kernel[j];
            }
        }
        result[i] = if val.is_finite() { val } else { 0.0 };
    }

    result
}

const AXIAL_SEARCH_MAX_STEPS: usize = 4096;

fn estimate_axial_length(signal: &[f64], center_idx: usize, axial_resolution: f64) -> f64 {
    if signal.is_empty() {
        return 0.0;
    }

    let safe_center = center_idx.min(signal.len().saturating_sub(1));
    let peak = signal[safe_center].abs();
    if peak < 1e-12 || !peak.is_finite() {
        return 0.0;
    }

    let threshold = peak * 0.5;

    let mut left = safe_center;
    let mut steps = 0;
    while left > 0 && signal[left].abs() > threshold && steps < AXIAL_SEARCH_MAX_STEPS {
        left -= 1;
        steps += 1;
    }

    let mut right = safe_center;
    steps = 0;
    while right < signal.len().saturating_sub(1)
        && signal[right].abs() > threshold
        && steps < AXIAL_SEARCH_MAX_STEPS
    {
        right += 1;
        steps += 1;
    }

    let width_samples = right.saturating_sub(left).min(signal.len());
    let width = width_samples as f64 * axial_resolution;
    if width.is_finite() { width } else { 0.0 }
}

pub fn compute_statistics(
    segments: &[PipeScanSegment],
    _wall_thickness: f64,
) -> crate::types::InversionResult {
    let mut all_depths: Vec<f64> = Vec::new();
    let mut total_critical = 0;
    let mut total_severe = 0;
    let mut total_moderate = 0;
    let mut total_mild = 0;

    let mut defect_map: Vec<Vec<DefectProfile>> = Vec::new();

    for segment in segments {
        for sensor_row in &segment.defect_map {
            for profile in sensor_row {
                let d = profile.depth_m;
                if d.is_finite() && d > 0.0 {
                    all_depths.push(d);
                }

                match profile.severity {
                    DefectSeverity::Critical => total_critical += 1,
                    DefectSeverity::Severe => total_severe += 1,
                    DefectSeverity::Moderate => total_moderate += 1,
                    DefectSeverity::Mild => total_mild += 1,
                    _ => {}
                }
            }
        }
    }

    if let Some(first) = segments.first() {
        let num_sensors = first.defect_map.len();
        for segment in segments {
            if segment.defect_map.is_empty() {
                continue;
            }
            let num_samples = segment.defect_map[0].len();

            if defect_map.is_empty() {
                defect_map = vec![Vec::with_capacity(all_depths.len() / num_sensors.max(1)); num_sensors];
            }

            for s in 0..num_sensors.min(segment.defect_map.len()) {
                defect_map[s].extend_from_slice(&segment.defect_map[s]);
            }
        }
    }

    let max_depth = all_depths
        .iter()
        .copied()
        .fold(0.0f64, f64::max);

    let avg_depth = if all_depths.is_empty() {
        0.0
    } else {
        all_depths.iter().sum::<f64>() / all_depths.len() as f64
    };

    let total_length = segments
        .last()
        .map(|s| s.end_position_m)
        .unwrap_or(0.0);

    let num_sensors = defect_map.len();
    let num_axial_points = if num_sensors > 0 { defect_map[0].len() } else { 0 };

    crate::types::InversionResult {
        total_length_m: total_length,
        num_sensors,
        num_axial_points,
        defect_map,
        max_depth_m: max_depth,
        avg_depth_m: avg_depth,
        critical_defect_count: total_critical,
        severe_defect_count: total_severe,
        moderate_defect_count: total_moderate,
        mild_defect_count: total_mild,
    }
}
