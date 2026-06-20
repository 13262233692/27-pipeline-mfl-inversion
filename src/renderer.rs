use std::fmt::Write as _;

use crate::types::{
    DefectProfile, DefectSeverity, InversionResult, WALL_THICKNESS_M, depth_to_ascii,
};

const DEFAULT_TERMINAL_WIDTH: usize = 120;
const DEFAULT_AXIAL_CHARS: usize = 100;
const HEADER_HEIGHT: usize = 3;
const LEGEND_WIDTH: usize = 30;

pub struct AsciiRenderer {
    terminal_width: usize,
    wall_thickness: f64,
    show_axis: bool,
    show_legend: bool,
    color_enabled: bool,
}

impl Default for AsciiRenderer {
    fn default() -> Self {
        Self {
            terminal_width: DEFAULT_TERMINAL_WIDTH,
            wall_thickness: WALL_THICKNESS_M,
            show_axis: true,
            show_legend: true,
            color_enabled: true,
        }
    }
}

impl AsciiRenderer {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_terminal_width(mut self, width: usize) -> Self {
        self.terminal_width = width;
        self
    }

    pub fn with_wall_thickness(mut self, thickness: f64) -> Self {
        self.wall_thickness = thickness;
        self
    }

    pub fn with_color(mut self, enabled: bool) -> Self {
        self.color_enabled = enabled;
        self
    }

    pub fn render_unfolded_map(&self, result: &InversionResult) -> String {
        let mut output = String::new();

        self.render_header(&mut output, result);
        self.render_defect_map(&mut output, result);
        self.render_legend(&mut output);
        self.render_summary(&mut output, result);

        output
    }

    fn render_header(&self, output: &mut String, result: &InversionResult) {
        let title = "╔══ MFL PIPELINE DEFECT INVERSION - UNFOLDED WALL MAP ══╗";
        let pad_len = (self.terminal_width.saturating_sub(title.len())) / 2;
        writeln!(output, "{:>pad$}", "", pad = pad_len + title.len()).unwrap();
        writeln!(output, "{:>pad$}", "", pad = pad_len + title.len()).unwrap();

        let info_line = format!(
            "  Length: {:.2}m | Sensors: {} | Axial Points: {} | Wall Thickness: {:.1}mm",
            result.total_length_m,
            result.num_sensors,
            result.num_axial_points,
            self.wall_thickness * 1000.0
        );

        let map_width = self.terminal_width.saturating_sub(LEGEND_WIDTH + 4);
        let _ = map_width;

        writeln!(output, "{}", info_line).unwrap();
        writeln!(output).unwrap();
    }

    fn render_defect_map(&self, output: &mut String, result: &InversionResult) {
        if result.num_sensors == 0 || result.num_axial_points == 0 {
            writeln!(output, "  [No data to display]").unwrap();
            return;
        }

        let axial_chars = (self.terminal_width - LEGEND_WIDTH - 8).min(DEFAULT_AXIAL_CHARS * 2);
        let axial_step = (result.num_axial_points / axial_chars.max(1)).max(1);
        let display_axial = (result.num_axial_points / axial_step).max(1);

        let sensor_step = (result.num_sensors / 16).max(1);
        let display_sensors = result.num_sensors / sensor_step;

        writeln!(
            output,
            "  ┌{:─<width$}┐  Defect Legend",
            "",
            width = display_axial
        ).unwrap();

        for sensor_row in 0..display_sensors {
            let sensor_idx = sensor_row * sensor_step;

            let row_label = format!("S{:02}", sensor_idx);
            write!(output, "  │").unwrap();

            for axial_idx in 0..display_axial {
                let actual_axial = axial_idx * axial_step;
                if actual_axial < result.defect_map[sensor_idx].len() {
                    let profile = &result.defect_map[sensor_idx][actual_axial];
                    let depth_ratio = profile.depth_m / self.wall_thickness;
                    let ch = depth_to_ascii(depth_ratio);

                    if self.color_enabled {
                        let color_code = profile.severity.color_code();
                        write!(output, "\x1b[38;5;{}m{}\x1b[0m", color_code, ch).unwrap();
                    } else {
                        write!(output, "{}", ch).unwrap();
                    }
                } else {
                    write!(output, " ").unwrap();
                }
            }

            let severity_label = self.get_row_severity_label(sensor_idx, sensor_step, result);
            writeln!(output, "│ {} {}", row_label, severity_label).unwrap();
        }

        writeln!(
            output,
            "  └{:─<width$}┘",
            "",
            width = display_axial
        ).unwrap();

        self.render_axial_ticks(output, result, display_axial, axial_step);
    }

    fn get_row_severity_label(
        &self,
        sensor_idx: usize,
        sensor_step: usize,
        result: &InversionResult,
    ) -> String {
        let mut max_severity = DefectSeverity::Normal;
        let mut max_depth = 0.0f64;

        for s in sensor_idx..(sensor_idx + sensor_step).min(result.num_sensors) {
            for profile in &result.defect_map[s] {
                if profile.depth_m > max_depth {
                    max_depth = profile.depth_m;
                    max_severity = profile.severity;
                }
            }
        }

        if self.color_enabled {
            format!(
                "\x1b[38;5;{}m{:8}\x1b[0m ({:.2}mm)",
                max_severity.color_code(),
                max_severity.label(),
                max_depth * 1000.0
            )
        } else {
            format!("{:8} ({:.2}mm)", max_severity.label(), max_depth * 1000.0)
        }
    }

    fn render_axial_ticks(
        &self,
        output: &mut String,
        result: &InversionResult,
        display_axial: usize,
        axial_step: usize,
    ) {
        let num_ticks = 5;
        let step = display_axial / num_ticks.max(1);

        write!(output, "   ").unwrap();

        for i in 0..=num_ticks {
            let pos = i * step;
            if pos < display_axial {
                let actual_sample = pos * axial_step;
                let position_m = actual_sample as f64 / result.num_axial_points as f64 * result.total_length_m;

                let label = format!("{:.0}m", position_m);
                if i > 0 {
                    let spaces = step.saturating_sub(label.len());
                    write!(output, "{:>width$}", label, width = label.len() + spaces).unwrap();
                } else {
                    write!(output, "{}", label).unwrap();
                }
            }
        }
        writeln!(output).unwrap();
        writeln!(output).unwrap();
    }

    fn render_legend(&self, output: &mut String) {
        writeln!(output, "  ═══════════════════ SEVERITY LEGEND ═══════════════════").unwrap();
        writeln!(output).unwrap();

        let severities = [
            (DefectSeverity::Normal, "0% - 10% wall loss"),
            (DefectSeverity::Mild, "10% - 25% wall loss"),
            (DefectSeverity::Moderate, "25% - 40% wall loss"),
            (DefectSeverity::Severe, "40% - 60% wall loss"),
            (DefectSeverity::Critical, "> 60% wall loss"),
        ];

        for (severity, desc) in &severities {
            let block = "████████";
            if self.color_enabled {
                writeln!(
                    output,
                    "    \x1b[38;5;{}m{}\x1b[0m  {:<8}  {}",
                    severity.color_code(),
                    block,
                    severity.label(),
                    desc
                ).unwrap();
            } else {
                writeln!(
                    output,
                    "    {}  {:<8}  {}",
                    block,
                    severity.label(),
                    desc
                ).unwrap();
            }
        }

        writeln!(output).unwrap();
        writeln!(output, "  Depth gradient:").unwrap();

        let gradient_chars = super::types::ASCII_GRADIENT;
        write!(output, "    ").unwrap();
        for (i, ch) in gradient_chars.iter().enumerate() {
            let ratio = i as f64 / (gradient_chars.len() - 1) as f64;
            let severity = DefectSeverity::from_depth_ratio(ratio);
            if self.color_enabled {
                write!(output, "\x1b[38;5;{}m{}\x1b[0m", severity.color_code(), ch).unwrap();
            } else {
                write!(output, "{}", ch).unwrap();
            }
        }
        writeln!(output, "  shallow → deep").unwrap();
        writeln!(output).unwrap();
    }

    fn render_summary(&self, output: &mut String, result: &InversionResult) {
        writeln!(output, "  ════════════════════ STATISTICS ══════════════════════").unwrap();
        writeln!(output).unwrap();

        writeln!(output, "    Total scan length     : {:.2} m", result.total_length_m).unwrap();
        writeln!(output, "    Number of sensors     : {}", result.num_sensors).unwrap();
        writeln!(output, "    Axial data points     : {}", result.num_axial_points).unwrap();
        writeln!(output, "    Maximum depth         : {:.3} mm ({:.1}% wall loss)",
            result.max_depth_m * 1000.0,
            result.max_depth_m / self.wall_thickness * 100.0
        ).unwrap();
        writeln!(output, "    Average depth         : {:.3} mm",
            result.avg_depth_m * 1000.0
        ).unwrap();
        writeln!(output).unwrap();

        let total_points = result.num_sensors * result.num_axial_points;
        let total_defects = result.critical_defect_count
            + result.severe_defect_count
            + result.moderate_defect_count
            + result.mild_defect_count;

        writeln!(output, "    Defect distribution by severity:").unwrap();
        writeln!(output).unwrap();

        if self.color_enabled {
            writeln!(
                output,
                "      \x1b[38;5;{}mCRITICAL\x1b[0m : {} points ({:.2}%)",
                DefectSeverity::Critical.color_code(),
                result.critical_defect_count,
                result.critical_defect_count as f64 / total_points.max(1) as f64 * 100.0
            ).unwrap();
            writeln!(
                output,
                "      \x1b[38;5;{}mSEVERE \x1b[0m : {} points ({:.2}%)",
                DefectSeverity::Severe.color_code(),
                result.severe_defect_count,
                result.severe_defect_count as f64 / total_points.max(1) as f64 * 100.0
            ).unwrap();
            writeln!(
                output,
                "      \x1b[38;5;{}mMODERATE\x1b[0m: {} points ({:.2}%)",
                DefectSeverity::Moderate.color_code(),
                result.moderate_defect_count,
                result.moderate_defect_count as f64 / total_points.max(1) as f64 * 100.0
            ).unwrap();
            writeln!(
                output,
                "      \x1b[38;5;{}mMILD   \x1b[0m : {} points ({:.2}%)",
                DefectSeverity::Mild.color_code(),
                result.mild_defect_count,
                result.mild_defect_count as f64 / total_points.max(1) as f64 * 100.0
            ).unwrap();
        } else {
            writeln!(
                output,
                "      CRITICAL : {} points ({:.2}%)",
                result.critical_defect_count,
                result.critical_defect_count as f64 / total_points.max(1) as f64 * 100.0
            ).unwrap();
            writeln!(
                output,
                "      SEVERE   : {} points ({:.2}%)",
                result.severe_defect_count,
                result.severe_defect_count as f64 / total_points.max(1) as f64 * 100.0
            ).unwrap();
            writeln!(
                output,
                "      MODERATE : {} points ({:.2}%)",
                result.moderate_defect_count,
                result.moderate_defect_count as f64 / total_points.max(1) as f64 * 100.0
            ).unwrap();
            writeln!(
                output,
                "      MILD     : {} points ({:.2}%)",
                result.mild_defect_count,
                result.mild_defect_count as f64 / total_points.max(1) as f64 * 100.0
            ).unwrap();
        }

        writeln!(output).unwrap();
        writeln!(
            output,
            "    Total anomalous points: {} ({:.2}% of scan)",
            total_defects,
            total_defects as f64 / total_points.max(1) as f64 * 100.0
        ).unwrap();
        writeln!(output).unwrap();
        writeln!(output, "  ══════════════════════════════════════════════════════").unwrap();
    }
}

pub fn render_defect_detail(profile: &DefectProfile, wall_thickness: f64, color: bool) -> String {
    let mut output = String::new();

    let depth_mm = profile.depth_m * 1000.0;
    let depth_pct = profile.depth_m / wall_thickness * 100.0;

    if color {
        writeln!(
            output,
            "\x1b[38;5;{}m[{}]\x1b[0m  Depth: {:.3}mm ({:.1}%)  Axial Length: {:.3}mm",
            profile.severity.color_code(),
            profile.severity.label(),
            depth_mm,
            depth_pct,
            profile.axial_length_m * 1000.0
        ).unwrap();
    } else {
        writeln!(
            output,
            "[{}]  Depth: {:.3}mm ({:.1}%)  Axial Length: {:.3}mm",
            profile.severity.label(),
            depth_mm,
            depth_pct,
            profile.axial_length_m * 1000.0
        ).unwrap();
    }

    output
}
