use crate::types::{DefectProfile, InversionResult};

const SAFETY_FACTOR: f64 = 1.5;

#[derive(Debug, Clone, Copy)]
pub struct AsmeB31gParams {
    pub smys_pa: f64,
    pub nominal_wall_thickness_m: f64,
    pub outer_diameter_m: f64,
    pub operating_pressure_pa: f64,
}

impl Default for AsmeB31gParams {
    fn default() -> Self {
        Self {
            smys_pa: 290_000_000.0,
            nominal_wall_thickness_m: 0.0127,
            outer_diameter_m: 0.3048,
            operating_pressure_pa: 10_000_000.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct CorrosionCluster {
    pub max_depth_m: f64,
    pub axial_length_m: f64,
    pub circumferential_width_m: f64,
    pub area_m2: f64,
    pub start_axial_m: f64,
    pub end_axial_m: f64,
    pub start_sensor: usize,
    pub end_sensor: usize,
}

#[derive(Debug, Clone)]
pub struct AsmeB31gResult {
    pub maop_pa: f64,
    pub burst_pressure_pa: f64,
    pub flow_stress_pa: f64,
    pub folias_factor_m: f64,
    pub corrosion_area_m2: f64,
    pub nominal_area_m2: f64,
    pub area_ratio: f64,
    pub max_corrosion: CorrosionCluster,
    pub is_overpressure: bool,
    pub safety_margin: f64,
}

pub fn evaluate_asme_b31g(
    result: &InversionResult,
    params: AsmeB31gParams,
) -> AsmeB31gResult {
    let max_cluster = find_max_corrosion_cluster(result, params.nominal_wall_thickness_m, params.outer_diameter_m);

    let l = max_cluster.axial_length_m.max(0.001);
    let d = max_cluster.max_depth_m;
    let t = params.nominal_wall_thickness_m;
    let d_o = params.outer_diameter_m;

    let area_ratio = if d > 0.0 && l > 0.0 {
        d / t
    } else {
        0.0
    };

    let sqrt_term = l * l / (d_o * t);
    let m = (1.0 + 0.8 * sqrt_term).sqrt();

    let flow_stress = params.smys_pa + 68_950_000.0;

    let one_minus_a_a0 = 1.0 - area_ratio;
    let one_minus_a_a0m = 1.0 - area_ratio / m;

    let burst_pressure = if one_minus_a_a0m.abs() < 1e-9 {
        0.0
    } else {
        2.0 * flow_stress * t * one_minus_a_a0 / (d_o * one_minus_a_a0m)
    };

    let maop = burst_pressure / SAFETY_FACTOR;

    let is_overpressure = params.operating_pressure_pa > maop;
    let safety_margin = maop / params.operating_pressure_pa;

    AsmeB31gResult {
        maop_pa: if maop.is_finite() { maop } else { 0.0 },
        burst_pressure_pa: if burst_pressure.is_finite() { burst_pressure } else { 0.0 },
        flow_stress_pa: flow_stress,
        folias_factor_m: m,
        corrosion_area_m2: max_cluster.area_m2,
        nominal_area_m2: t * l,
        area_ratio,
        max_corrosion: max_cluster,
        is_overpressure,
        safety_margin: if safety_margin.is_finite() { safety_margin } else { 0.0 },
    }
}

fn find_max_corrosion_cluster(
    result: &InversionResult,
    wall_thickness: f64,
    outer_diameter: f64,
) -> CorrosionCluster {
    let num_sensors = result.defect_map.len();
    let num_axial = if num_sensors > 0 { result.defect_map[0].len() } else { 0 };

    if num_sensors == 0 || num_axial == 0 {
        return CorrosionCluster {
            max_depth_m: 0.0,
            axial_length_m: 0.0,
            circumferential_width_m: 0.0,
            area_m2: 0.0,
            start_axial_m: 0.0,
            end_axial_m: 0.0,
            start_sensor: 0,
            end_sensor: 0,
        };
    }

    let threshold_depth = wall_thickness * 0.1;
    let axial_resolution = result.total_length_m / num_axial as f64;
    let circumferential_resolution = std::f64::consts::PI * outer_diameter / num_sensors as f64;

    let mut visited = vec![vec![false; num_axial]; num_sensors];
    let mut max_area = 0.0f64;
    let mut max_cluster = CorrosionCluster {
        max_depth_m: 0.0,
        axial_length_m: 0.0,
        circumferential_width_m: 0.0,
        area_m2: 0.0,
        start_axial_m: 0.0,
        end_axial_m: 0.0,
        start_sensor: 0,
        end_sensor: 0,
    };

    for sensor in 0..num_sensors {
        for axial in 0..num_axial {
            if visited[sensor][axial] || result.defect_map[sensor][axial].depth_m <= threshold_depth {
                continue;
            }

            let (cluster, area) = flood_fill_cluster(
                &result.defect_map,
                sensor,
                axial,
                &mut visited,
                threshold_depth,
                axial_resolution,
                circumferential_resolution,
                wall_thickness,
            );

            if area > max_area {
                max_area = area;
                max_cluster = cluster;
            }
        }
    }

    max_cluster
}

fn flood_fill_cluster(
    defect_map: &[Vec<DefectProfile>],
    start_sensor: usize,
    start_axial: usize,
    visited: &mut [Vec<bool>],
    threshold: f64,
    axial_res: f64,
    circ_res: f64,
    _wall_thickness: f64,
) -> (CorrosionCluster, f64) {
    let num_sensors = defect_map.len();
    let num_axial = defect_map[0].len();

    let mut stack = vec![(start_sensor, start_axial)];
    let mut min_sensor = start_sensor;
    let mut max_sensor = start_sensor;
    let mut min_axial = start_axial;
    let mut max_axial = start_axial;
    let mut total_depth = 0.0f64;
    let mut max_depth = 0.0f64;
    let mut count = 0usize;

    visited[start_sensor][start_axial] = true;

    while let Some((s, a)) = stack.pop() {
        let d = defect_map[s][a].depth_m;
        if d <= threshold {
            continue;
        }

        total_depth += d;
        count += 1;
        max_depth = max_depth.max(d);

        min_sensor = min_sensor.min(s);
        max_sensor = max_sensor.max(s);
        min_axial = min_axial.min(a);
        max_axial = max_axial.max(a);

        let offsets = [
            (s as isize - 1, a as isize),
            (s as isize + 1, a as isize),
            (s as isize, a as isize - 1),
            (s as isize, a as isize + 1),
            (s as isize - 1, a as isize - 1),
            (s as isize + 1, a as isize - 1),
            (s as isize - 1, a as isize + 1),
            (s as isize + 1, a as isize + 1),
        ];

        for &(ds, da) in &offsets {
            if ds < 0 || da < 0 {
                continue;
            }
            let ns = ds as usize;
            let na = da as usize;
            if ns < num_sensors && na < num_axial && !visited[ns][na] && defect_map[ns][na].depth_m > threshold {
                visited[ns][na] = true;
                stack.push((ns, na));
            }
        }
    }

    let avg_depth = if count > 0 { total_depth / count as f64 } else { 0.0 };
    let axial_length = (max_axial - min_axial + 1) as f64 * axial_res;
    let circ_width = (max_sensor - min_sensor + 1) as f64 * circ_res;
    let area = axial_length * avg_depth;

    let cluster = CorrosionCluster {
        max_depth_m: max_depth,
        axial_length_m: axial_length,
        circumferential_width_m: circ_width,
        area_m2: area,
        start_axial_m: min_axial as f64 * axial_res,
        end_axial_m: (max_axial + 1) as f64 * axial_res,
        start_sensor: min_sensor,
        end_sensor: max_sensor,
    };

    (cluster, area)
}

pub fn render_b31g_report(
    result: &AsmeB31gResult,
    params: &AsmeB31gParams,
    color: bool,
) -> String {
    let mut out = String::new();

    if color {
        out.push_str("\n  \x1b[38;5;39m╔═══════════════════════════════════════════════════════════════════╗\x1b[0m\n");
        out.push_str("  \x1b[38;5;39m║\x1b[0m       \x1b[1;38;5;208mASME B31G  REMAINING STRENGTH ASSESSMENT  \x1b[0m           \x1b[38;5;39m║\x1b[0m\n");
        out.push_str("  \x1b[38;5;39m╚═══════════════════════════════════════════════════════════════════╝\x1b[0m\n");
    } else {
        out.push_str("\n  ════════════════════════════════════════════════════════════════\n");
        out.push_str("        ASME B31G  REMAINING STRENGTH ASSESSMENT\n");
        out.push_str("  ════════════════════════════════════════════════════════════════\n");
    }

    out.push_str(&format!("  Pipe OD:           {:.1} mm\n", params.outer_diameter_m * 1000.0));
    out.push_str(&format!("  Nominal wall:      {:.1} mm\n", params.nominal_wall_thickness_m * 1000.0));
    out.push_str(&format!("  SMYS:              {:.0} MPa\n", params.smys_pa / 1e6));
    out.push_str(&format!("  Operating press:   {:.2} MPa\n", params.operating_pressure_pa / 1e6));
    out.push('\n');

    out.push_str("  ───────────────────────────────────────────────────────────────\n");
    out.push_str("  MAX CORROSION CLUSTER\n");
    out.push_str("  ───────────────────────────────────────────────────────────────\n");
    out.push_str(&format!("  Max depth:         {:.2} mm ({:.1}% wall loss)\n",
        result.max_corrosion.max_depth_m * 1000.0,
        result.max_corrosion.max_depth_m / params.nominal_wall_thickness_m * 100.0));
    out.push_str(&format!("  Axial length:      {:.2} m\n", result.max_corrosion.axial_length_m));
    out.push_str(&format!("  Circ width:        {:.2} m\n", result.max_corrosion.circumferential_width_m));
    out.push_str(&format!("  Location:          {:.2} m - {:.2} m\n",
        result.max_corrosion.start_axial_m,
        result.max_corrosion.end_axial_m));
    out.push_str(&format!("  Sensors:           S{:02} - S{:02}\n",
        result.max_corrosion.start_sensor,
        result.max_corrosion.end_sensor));
    out.push_str(&format!("  Area ratio (A/A0): {:.3}\n", result.area_ratio));
    out.push('\n');

    out.push_str("  ───────────────────────────────────────────────────────────────\n");
    out.push_str("  BURST PRESSURE CALCULATION\n");
    out.push_str("  ───────────────────────────────────────────────────────────────\n");
    out.push_str(&format!("  Flow stress:       {:.0} MPa\n", result.flow_stress_pa / 1e6));
    out.push_str(&format!("  Folias factor M:   {:.3}\n", result.folias_factor_m));
    out.push_str(&format!("  Burst pressure:    {:.2} MPa\n", result.burst_pressure_pa / 1e6));
    out.push_str(&format!("  MAOP:              {:.2} MPa\n", result.maop_pa / 1e6));
    out.push_str(&format!("  Safety factor FOS: {:.2}\n", SAFETY_FACTOR));
    out.push('\n');

    out.push_str("  ───────────────────────────────────────────────────────────────\n");
    out.push_str("  SAFETY ASSESSMENT\n");
    out.push_str("  ───────────────────────────────────────────────────────────────\n");
    out.push_str(&format!("  Operating press:   {:.2} MPa\n", params.operating_pressure_pa / 1e6));
    out.push_str(&format!("  Allowable MAOP:    {:.2} MPa\n", result.maop_pa / 1e6));
    out.push_str(&format!("  Safety margin:     {:.2}x\n", result.safety_margin));

    out
}

pub fn render_critical_alert() -> String {
    let mut s = String::new();
    s.push_str("\n\n");
    s.push_str("  \x1b[41;1;37m╔══════════════════════════════════════════════════════════════════════════════════╗\x1b[0m\n");
    s.push_str("  \x1b[41;1;37m║\x1b[0m  \x1b[41;1;37m▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\x1b[0m\n");
    s.push_str("  \x1b[41;1;37m║\x1b[0m  \x1b[41;1;37m▓▓  【 极 高 穿 孔 爆 破 风 险 】  EXTREME PERFORATION & BURST HAZARD  ▓▓\x1b[0m\n");
    s.push_str("  \x1b[41;1;37m║\x1b[0m  \x1b[41;1;37m▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓\x1b[0m\n");
    s.push_str("  \x1b[41;1;37m║\x1b[0m                                                                              \x1b[41;1;37m║\x1b[0m\n");
    s.push_str("  \x1b[41;1;37m║\x1b[0m  \x1b[41;1;37m  OPERATING PRESSURE EXCEEDS MAOP LIMIT - IMMEDIATE SHUTDOWN REQUIRED  \x1b[0m\n");
    s.push_str("  \x1b[41;1;37m║\x1b[0m  \x1b[41;1;37m  运 行 压 力 已 超 过 MAOP 临 界 线  -  立 即 停 运 隔 离 管 段  \x1b[0m\n");
    s.push_str("  \x1b[41;1;37m║\x1b[0m                                                                              \x1b[41;1;37m║\x1b[0m\n");
    s.push_str("  \x1b[41;1;37m╚══════════════════════════════════════════════════════════════════════════════════╝\x1b[0m\n");
    s.push_str("\n\n");
    s
}
