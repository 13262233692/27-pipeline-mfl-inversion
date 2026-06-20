use std::fs::File;
use std::path::Path;

use anyhow::{Context, Result};
use byteorder::{ByteOrder as ByteOrderTrait, LittleEndian, BigEndian};
use memmap2::Mmap;
use rayon::prelude::*;

use crate::types::{FileFormat, MflSample, PipeScanSegment, ByteOrder};
use crate::dipole::ProgressReporter;

const CHUNK_SAMPLES: usize = 10_000;
const OVERLAP_SAMPLES: usize = 100;
const READ_PROGRESS_THROTTLE: u64 = 128;

pub struct MflFileReader {
    format: FileFormat,
    mmap: Mmap,
    file_size: usize,
    total_samples: usize,
    num_sensors: usize,
}

impl MflFileReader {
    pub fn open<P: AsRef<Path>>(path: P, format: FileFormat) -> Result<Self> {
        let file = File::open(path.as_ref())
            .with_context(|| format!("Failed to open MFL file: {}", path.as_ref().display()))?;
        let file_size = file.metadata()?.len() as usize;
        let mmap = unsafe { Mmap::map(&file) }
            .context("Failed to memory-map MFL file")?;

        let bytes_per_frame = format.num_sensors * 3 * format.bytes_per_sample;
        let total_samples = file_size / bytes_per_frame;

        Ok(Self {
            format,
            mmap,
            file_size,
            total_samples,
            num_sensors: format.num_sensors,
        })
    }

    pub fn total_samples(&self) -> usize {
        self.total_samples
    }

    pub fn num_sensors(&self) -> usize {
        self.num_sensors
    }

    pub fn file_size(&self) -> usize {
        self.file_size
    }

    pub fn total_length_m(&self) -> f64 {
        self.total_samples as f64 / self.format.sample_rate_hz * self.format.pig_speed_m_s
    }

    pub fn read_all_parallel(&self) -> Result<Vec<PipeScanSegment>> {
        self.read_all_parallel_with_progress(None)
    }

    pub fn read_all_parallel_with_progress(&self, reporter: Option<&ProgressReporter>) -> Result<Vec<PipeScanSegment>> {
        let chunk_size = CHUNK_SAMPLES;
        let num_chunks = (self.total_samples + chunk_size - 1) / chunk_size;

        let segments: Result<Vec<PipeScanSegment>> = (0..num_chunks)
            .into_par_iter()
            .map(|chunk_idx| {
                let start_sample = chunk_idx * chunk_size;
                let end_sample = ((chunk_idx + 1) * chunk_size + OVERLAP_SAMPLES).min(self.total_samples);
                let actual_start = start_sample.saturating_sub(if chunk_idx > 0 { OVERLAP_SAMPLES } else { 0 });
                let seg = self.read_segment(actual_start, end_sample)?;

                if let Some(r) = reporter {
                    let work_units = ((end_sample - actual_start) * self.num_sensors) as u64;
                    let mut tick_accum = work_units;
                    while tick_accum >= READ_PROGRESS_THROTTLE {
                        r.tick_read(READ_PROGRESS_THROTTLE);
                        tick_accum -= READ_PROGRESS_THROTTLE;
                    }
                    if tick_accum > 0 {
                        r.tick_read(tick_accum);
                    }
                }

                Ok(seg)
            })
            .collect();

        segments
    }

    fn read_segment(&self, start_sample: usize, end_sample: usize) -> Result<PipeScanSegment> {
        let num_samples = end_sample - start_sample;
        let bytes_per_sample = self.format.bytes_per_sample;
        let bytes_per_frame = self.num_sensors * 3 * bytes_per_sample;

        let start_byte = start_sample * bytes_per_frame;
        let end_byte = end_sample * bytes_per_frame;

        if end_byte > self.file_size {
            anyhow::bail!(
                "Segment out of bounds: end_byte={}, file_size={}",
                end_byte, self.file_size
            );
        }

        let data = &self.mmap[start_byte..end_byte];
        let mut segment = PipeScanSegment::new(
            self.sample_to_position(start_sample),
            self.sample_to_position(end_sample),
            self.num_sensors,
            num_samples,
        );

        match self.format.byte_order {
            ByteOrder::LittleEndian => {
                self.decode_samples_le(data, &mut segment.sensor_data, num_samples)
            }
            ByteOrder::BigEndian => {
                self.decode_samples_be(data, &mut segment.sensor_data, num_samples)
            }
        }

        Ok(segment)
    }

    fn decode_samples_le(
        &self,
        data: &[u8],
        sensor_data: &mut [Vec<MflSample>],
        num_samples: usize,
    ) {
        let bps = self.format.bytes_per_sample;
        let frame_size = self.num_sensors * 3 * bps;

        for s in 0..num_samples {
            let frame_offset = s * frame_size;
            for sensor in 0..self.num_sensors {
                let sensor_offset = frame_offset + sensor * 3 * bps;

                let axial = match bps {
                    2 => LittleEndian::read_i16(&data[sensor_offset..sensor_offset + 2]) as f64,
                    _ => LittleEndian::read_i32(&data[sensor_offset..sensor_offset + 4]) as f64,
                };

                let radial_offset = sensor_offset + bps;
                let radial = match bps {
                    2 => LittleEndian::read_i16(&data[radial_offset..radial_offset + 2]) as f64,
                    _ => LittleEndian::read_i32(&data[radial_offset..radial_offset + 4]) as f64,
                };

                let circ_offset = sensor_offset + 2 * bps;
                let circumferential = match bps {
                    2 => LittleEndian::read_i16(&data[circ_offset..circ_offset + 2]) as f64,
                    _ => LittleEndian::read_i32(&data[circ_offset..circ_offset + 4]) as f64,
                };

                sensor_data[sensor][s] = MflSample {
                    axial: axial / 1e6 / COIL_SENSITIVITY_V_T,
                    radial: radial / 1e6 / COIL_SENSITIVITY_V_T,
                    circumferential: circumferential / 1e6 / COIL_SENSITIVITY_V_T,
                };
            }
        }
    }

    fn decode_samples_be(
        &self,
        data: &[u8],
        sensor_data: &mut [Vec<MflSample>],
        num_samples: usize,
    ) {
        let bps = self.format.bytes_per_sample;
        let frame_size = self.num_sensors * 3 * bps;

        for s in 0..num_samples {
            let frame_offset = s * frame_size;
            for sensor in 0..self.num_sensors {
                let sensor_offset = frame_offset + sensor * 3 * bps;

                let axial = match bps {
                    2 => BigEndian::read_i16(&data[sensor_offset..sensor_offset + 2]) as f64,
                    _ => BigEndian::read_i32(&data[sensor_offset..sensor_offset + 4]) as f64,
                };

                let radial_offset = sensor_offset + bps;
                let radial = match bps {
                    2 => BigEndian::read_i16(&data[radial_offset..radial_offset + 2]) as f64,
                    _ => BigEndian::read_i32(&data[radial_offset..radial_offset + 4]) as f64,
                };

                let circ_offset = sensor_offset + 2 * bps;
                let circumferential = match bps {
                    2 => BigEndian::read_i16(&data[circ_offset..circ_offset + 2]) as f64,
                    _ => BigEndian::read_i32(&data[circ_offset..circ_offset + 4]) as f64,
                };

                sensor_data[sensor][s] = MflSample {
                    axial: axial / 1e6 / COIL_SENSITIVITY_V_T,
                    radial: radial / 1e6 / COIL_SENSITIVITY_V_T,
                    circumferential: circumferential / 1e6 / COIL_SENSITIVITY_V_T,
                };
            }
        }
    }

    fn sample_to_position(&self, sample_idx: usize) -> f64 {
        sample_idx as f64 / self.format.sample_rate_hz * self.format.pig_speed_m_s
    }

    pub fn stream_chunks<F>(&self, mut callback: F) -> Result<()>
    where
        F: FnMut(PipeScanSegment) -> Result<()>,
    {
        let segments = self.read_all_parallel()?;
        for seg in segments {
            callback(seg)?;
        }
        Ok(())
    }
}

pub fn generate_test_file<P: AsRef<Path>>(
    path: P,
    num_samples: usize,
    num_sensors: usize,
) -> Result<()> {
    use std::io::Write;

    let mut file = File::create(path.as_ref())?;

    let pipe_radius = 0.1524;

    for s in 0..num_samples {
        let axial_pos = s as f64 / SAMPLE_RATE_HZ * PIG_SPEED_M_S;

        for sensor in 0..num_sensors {
            let sensor_angle = 2.0 * std::f64::consts::PI * sensor as f64 / num_sensors as f64;

            let (defect_axial, defect_depth) = if s >= 2000 && s <= 2300 && sensor == 32 {
                let defect_center = 2150;
                let dist_from_center = (s as isize - defect_center).abs() as f64;
                let gaussian = (-dist_from_center * dist_from_center / (50.0 * 50.0)).exp();
                (500e-6 * gaussian, 0.005 * gaussian)
            } else if s >= 5000 && s <= 5600 && sensor >= 10 && sensor <= 20 {
                let defect_center = 5300;
                let dist_from_center = (s as isize - defect_center).abs() as f64;
                let gaussian = (-dist_from_center * dist_from_center / (150.0 * 150.0)).exp();
                let sensor_mid = 15.0;
                let sensor_dist = (sensor as f64 - sensor_mid).abs();
                let sensor_gauss = (-sensor_dist * sensor_dist / (3.0 * 3.0)).exp();
                (800e-6 * gaussian * sensor_gauss, 0.008 * gaussian * sensor_gauss)
            } else if s >= 8000 && s <= 8500 && sensor >= 45 && sensor <= 55 {
                let defect_center = 8250;
                let dist_from_center = (s as isize - defect_center).abs() as f64;
                let gaussian = (-dist_from_center * dist_from_center / (100.0 * 100.0)).exp();
                let sensor_mid = 50.0;
                let sensor_dist = (sensor as f64 - sensor_mid).abs();
                let sensor_gauss = (-sensor_dist * sensor_dist / (4.0 * 4.0)).exp();
                (1200e-6 * gaussian * sensor_gauss, 0.012 * gaussian * sensor_gauss)
            } else {
                (0.0, 0.0)
            };

            let noise_amp = 1e-6;
            let axial = defect_axial * 5.0 + noise_amp * fast_rand(s * 1000 + sensor) + 5e-6 * (axial_pos * 2.0).sin();
            let radial = -defect_depth * 5e-3 + noise_amp * fast_rand(s * 1000 + sensor + 1) + 2e-6 * (axial_pos * 3.5).cos();
            let circumferential = 0.5 * defect_axial * 5.0 * sensor_angle.sin()
                + noise_amp * fast_rand(s * 1000 + sensor + 2);

            let axial_volts = axial * COIL_SENSITIVITY_V_T;
            let radial_volts = radial * COIL_SENSITIVITY_V_T;
            let circumferential_volts = circumferential * COIL_SENSITIVITY_V_T;

            let axial_microvolts = (axial_volts * 1e6) as i32;
            let radial_microvolts = (radial_volts * 1e6) as i32;
            let circ_microvolts = (circumferential_volts * 1e6) as i32;

            let mut buf = [0u8; 12];
            LittleEndian::write_i32(&mut buf[0..4], axial_microvolts);
            LittleEndian::write_i32(&mut buf[4..8], radial_microvolts);
            LittleEndian::write_i32(&mut buf[8..12], circ_microvolts);
            file.write_all(&buf)?;
        }
    }

    Ok(())
}

fn fast_rand(seed: usize) -> f64 {
    let x = seed as u64;
    let x = x.wrapping_mul(6364136223846793005);
    let x = x.wrapping_add(1442695040888963407);
    let x = x ^ (x >> 21);
    (x as f64) / (u64::MAX as f64) * 2.0 - 1.0
}

use crate::types::{SAMPLE_RATE_HZ, PIG_SPEED_M_S, COIL_SENSITIVITY_V_T};
