use std::fmt;

pub const SENSOR_COUNT: usize = 64;
pub const SAMPLE_RATE_HZ: f64 = 1000.0;
pub const PIG_SPEED_M_S: f64 = 2.0;
pub const PIPE_DIAMETER_M: f64 = 0.3048;
pub const WALL_THICKNESS_M: f64 = 0.0127;
pub const MAGNETIC_PERMEABILITY: f64 = 4.0 * std::f64::consts::PI * 1e-7;
pub const COIL_SENSITIVITY_V_T: f64 = 100.0;

#[derive(Debug, Clone, Copy)]
pub struct MflSample {
    pub axial: f64,
    pub radial: f64,
    pub circumferential: f64,
}

impl MflSample {
    pub fn magnitude(&self) -> f64 {
        (self.axial * self.axial + self.radial * self.radial + self.circumferential * self.circumferential).sqrt()
    }
}

impl Default for MflSample {
    fn default() -> Self {
        Self {
            axial: 0.0,
            radial: 0.0,
            circumferential: 0.0,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct DefectProfile {
    pub depth_m: f64,
    pub axial_length_m: f64,
    pub severity: DefectSeverity,
}

impl Default for DefectProfile {
    fn default() -> Self {
        Self {
            depth_m: 0.0,
            axial_length_m: 0.0,
            severity: DefectSeverity::Normal,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub enum DefectSeverity {
    Normal = 0,
    Mild = 1,
    Moderate = 2,
    Severe = 3,
    Critical = 4,
}

impl DefectSeverity {
    pub fn from_depth_ratio(ratio: f64) -> Self {
        if ratio < 0.1 {
            DefectSeverity::Normal
        } else if ratio < 0.25 {
            DefectSeverity::Mild
        } else if ratio < 0.4 {
            DefectSeverity::Moderate
        } else if ratio < 0.6 {
            DefectSeverity::Severe
        } else {
            DefectSeverity::Critical
        }
    }

    pub fn color_code(&self) -> u8 {
        match self {
            DefectSeverity::Normal => 82,
            DefectSeverity::Mild => 190,
            DefectSeverity::Moderate => 214,
            DefectSeverity::Severe => 202,
            DefectSeverity::Critical => 196,
        }
    }

    pub fn label(&self) -> &'static str {
        match self {
            DefectSeverity::Normal => "NORMAL",
            DefectSeverity::Mild => "MILD",
            DefectSeverity::Moderate => "MODERATE",
            DefectSeverity::Severe => "SEVERE",
            DefectSeverity::Critical => "CRITICAL",
        }
    }
}

impl fmt::Display for DefectSeverity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.label())
    }
}

#[derive(Debug, Clone)]
pub struct PipeScanSegment {
    pub start_position_m: f64,
    pub end_position_m: f64,
    pub sensor_data: Vec<Vec<MflSample>>,
    pub defect_map: Vec<Vec<DefectProfile>>,
}

impl PipeScanSegment {
    pub fn new(start_pos: f64, end_pos: f64, num_sensors: usize, num_samples: usize) -> Self {
        let sensor_data = vec![vec![MflSample::default(); num_samples]; num_sensors];
        let defect_map = vec![vec![DefectProfile::default(); num_samples]; num_sensors];
        Self {
            start_position_m: start_pos,
            end_position_m: end_pos,
            sensor_data,
            defect_map,
        }
    }

    pub fn axial_resolution_m(&self) -> f64 {
        let num_samples = self.sensor_data[0].len();
        (self.end_position_m - self.start_position_m) / num_samples.max(1) as f64
    }
}

#[derive(Debug, Clone)]
pub struct InversionResult {
    pub total_length_m: f64,
    pub num_sensors: usize,
    pub num_axial_points: usize,
    pub defect_map: Vec<Vec<DefectProfile>>,
    pub max_depth_m: f64,
    pub avg_depth_m: f64,
    pub critical_defect_count: usize,
    pub severe_defect_count: usize,
    pub moderate_defect_count: usize,
    pub mild_defect_count: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct FileFormat {
    pub bytes_per_sample: usize,
    pub num_sensors: usize,
    pub sample_rate_hz: f64,
    pub pig_speed_m_s: f64,
    pub byte_order: ByteOrder,
}

#[derive(Debug, Clone, Copy)]
pub enum ByteOrder {
    LittleEndian,
    BigEndian,
}

impl Default for FileFormat {
    fn default() -> Self {
        Self {
            bytes_per_sample: 4,
            num_sensors: SENSOR_COUNT,
            sample_rate_hz: SAMPLE_RATE_HZ,
            pig_speed_m_s: PIG_SPEED_M_S,
            byte_order: ByteOrder::LittleEndian,
        }
    }
}

pub const ASCII_GRADIENT: [char; 10] = [' ', '·', '░', '▒', '▓', '█', '▇', '▆', '▅', '▃'];

pub fn depth_to_ascii(depth_ratio: f64) -> char {
    let idx = ((depth_ratio.clamp(0.0, 1.0)) * (ASCII_GRADIENT.len() - 1) as f64).round() as usize;
    ASCII_GRADIENT[idx.min(ASCII_GRADIENT.len() - 1)]
}
