use std::path::PathBuf;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use num_cpus;
use rayon;

mod types;
mod file_reader;
mod dipole;
mod renderer;

use crate::types::{FileFormat, ByteOrder, SENSOR_COUNT, SAMPLE_RATE_HZ, PIG_SPEED_M_S, WALL_THICKNESS_M};
use crate::file_reader::{MflFileReader, generate_test_file};
use crate::dipole::{DipoleInverter, compute_statistics};
use crate::renderer::AsciiRenderer;

#[derive(Parser, Debug)]
#[command(
    name = "mfl-inv",
    version = "0.1.0",
    about = "MFL (Magnetic Flux Leakage) Pipeline Defect Inversion CLI",
    long_about = "Hardcore CLI tool for offshore oil & gas pipeline inspection. \
                  Performs magnetic dipole inversion on raw MFL binary dump files \
                  from intelligent PIGs and outputs ANSI-colored ASCII wall defect maps.",
    author = "Pipeline Inspection Team"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    #[arg(short, long, default_value_t = false, help = "Disable ANSI color output")]
    no_color: bool,

    #[arg(short = 'j', long, help = "Number of parallel threads (default: all CPUs)")]
    threads: Option<usize>,
}

#[derive(Subcommand, Debug)]
enum Commands {
    #[command(name = "analyze", about = "Analyze an MFL data file and display defect map")]
    Analyze {
        #[arg(help = "Path to the MFL binary dump file")]
        input: PathBuf,

        #[arg(short = 's', long, default_value_t = SENSOR_COUNT, help = "Number of sensor channels")]
        sensors: usize,

        #[arg(short = 'r', long, default_value_t = SAMPLE_RATE_HZ, help = "Sample rate in Hz")]
        sample_rate: f64,

        #[arg(short = 'v', long, default_value_t = PIG_SPEED_M_S, help = "PIG speed in m/s")]
        pig_speed: f64,

        #[arg(short = 'w', long, default_value_t = WALL_THICKNESS_M, help = "Pipe wall thickness in meters")]
        wall_thickness: f64,

        #[arg(long, default_value = "little", help = "Byte order: little or big")]
        byte_order: String,

        #[arg(short = 'b', long, default_value_t = 4, help = "Bytes per sample (2, 4)")]
        bytes_per_sample: usize,

        #[arg(short = 'o', long, help = "Output file for JSON report")]
        output: Option<PathBuf>,
    },

    #[command(name = "gentest", about = "Generate a test MFL data file with simulated defects")]
    GenTest {
        #[arg(help = "Output path for the test file")]
        output: PathBuf,

        #[arg(short = 'n', long, default_value_t = 10000, help = "Number of samples to generate")]
        samples: usize,

        #[arg(short = 's', long, default_value_t = SENSOR_COUNT, help = "Number of sensor channels")]
        sensors: usize,
    },

    #[command(name = "info", about = "Display information about an MFL data file")]
    Info {
        #[arg(help = "Path to the MFL binary dump file")]
        input: PathBuf,

        #[arg(short = 's', long, default_value_t = SENSOR_COUNT, help = "Number of sensor channels")]
        sensors: usize,

        #[arg(short = 'b', long, default_value_t = 4, help = "Bytes per sample (2, 4)")]
        bytes_per_sample: usize,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    let thread_count = cli.threads.unwrap_or_else(num_cpus::get);
    rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build_global()
        .ok();

    let color_enabled = !cli.no_color && atty_stdout();

    match &cli.command {
        Commands::Analyze {
            input,
            sensors,
            sample_rate,
            pig_speed,
            wall_thickness,
            byte_order,
            bytes_per_sample,
            output,
        } => {
            run_analyze(
                input,
                *sensors,
                *sample_rate,
                *pig_speed,
                *wall_thickness,
                byte_order,
                *bytes_per_sample,
                output.as_deref(),
                color_enabled,
            )
        }
        Commands::GenTest {
            output,
            samples,
            sensors,
        } => run_gentest(output, *samples, *sensors, color_enabled),
        Commands::Info {
            input,
            sensors,
            bytes_per_sample,
        } => run_info(input, *sensors, *bytes_per_sample, color_enabled),
    }
}

fn run_analyze(
    input: &std::path::Path,
    num_sensors: usize,
    sample_rate: f64,
    pig_speed: f64,
    wall_thickness: f64,
    byte_order_str: &str,
    bytes_per_sample: usize,
    _output: Option<&std::path::Path>,
    color_enabled: bool,
) -> Result<()> {
    print_banner(color_enabled);

    let byte_order = match byte_order_str.to_lowercase().as_str() {
        "little" | "le" => ByteOrder::LittleEndian,
        "big" | "be" => ByteOrder::BigEndian,
        _ => anyhow::bail!("Invalid byte order: {}. Use 'little' or 'big'", byte_order_str),
    };

    let format = FileFormat {
        bytes_per_sample,
        num_sensors,
        sample_rate_hz: sample_rate,
        pig_speed_m_s: pig_speed,
        byte_order,
    };

    println!("  [CONFIG] Loading file: {}", input.display());
    println!("  [CONFIG] Sensors: {}, Sample rate: {} Hz", num_sensors, sample_rate);
    println!("  [CONFIG] PIG speed: {} m/s, Wall thickness: {:.3} mm", pig_speed, wall_thickness * 1000.0);
    println!("  [CONFIG] Byte order: {:?}, Bytes/sample: {}", byte_order, bytes_per_sample);
    println!();

    let start = Instant::now();
    let reader = MflFileReader::open(input, format)
        .with_context(|| format!("Failed to open MFL file: {}", input.display()))?;

    let file_size_mb = reader.file_size() as f64 / (1024.0 * 1024.0);
    println!("  [LOAD] File size: {:.2} MB", file_size_mb);
    println!("  [LOAD] Total samples: {}", reader.total_samples());
    println!("  [LOAD] Pipeline length: {:.2} m", reader.total_length_m());
    println!();

    println!("  [READ] Reading file in parallel chunks...");
    let mut segments = reader.read_all_parallel()?;
    let read_time = start.elapsed();
    println!("  [READ] Read complete in {:.2?}", read_time);
    println!();

    println!("  [INVERT] Running magnetic dipole inversion...");
    let inverter = DipoleInverter::new()
        .with_wall_thickness(wall_thickness);

    let invert_start = Instant::now();
    for segment in &mut segments {
        inverter.invert_segment_parallel(segment);
    }
    let invert_time = invert_start.elapsed();
    println!("  [INVERT] Inversion complete in {:.2?}", invert_time);
    println!();

    println!("  [STATS] Computing defect statistics...");
    let result = compute_statistics(&segments, wall_thickness);
    let stats_time = start.elapsed();
    println!("  [STATS] Statistics computed in {:.2?}", stats_time);
    println!();

    println!("  [RENDER] Generating ASCII wall defect map...");
    let renderer = AsciiRenderer::new()
        .with_wall_thickness(wall_thickness)
        .with_color(color_enabled);

    let render = renderer.render_unfolded_map(&result);
    println!();
    println!("{}", render);

    let total_time = start.elapsed();
    println!();
    println!("  [DONE] Total processing time: {:.2?}", total_time);
    println!();

    Ok(())
}

fn run_gentest(
    output: &std::path::Path,
    num_samples: usize,
    num_sensors: usize,
    color_enabled: bool,
) -> Result<()> {
    print_banner(color_enabled);

    println!("  [GEN] Generating test MFL data file...");
    println!("  [GEN] Output: {}", output.display());
    println!("  [GEN] Samples: {}, Sensors: {}", num_samples, num_sensors);

    let start = Instant::now();
    generate_test_file(output, num_samples, num_sensors)
        .with_context(|| format!("Failed to generate test file: {}", output.display()))?;

    let file_size = std::fs::metadata(output)?.len();
    let elapsed = start.elapsed();

    println!("  [GEN] File size: {:.2} MB", file_size as f64 / (1024.0 * 1024.0));
    println!("  [GEN] Generated in {:.2?}", elapsed);
    println!();
    println!("  Test file generated successfully!");
    println!("  Run 'mfl-inv analyze {}' to inspect it.", output.display());
    println!();

    Ok(())
}

fn run_info(
    input: &std::path::Path,
    num_sensors: usize,
    bytes_per_sample: usize,
    color_enabled: bool,
) -> Result<()> {
    print_banner(color_enabled);

    let metadata = std::fs::metadata(input)
        .with_context(|| format!("Failed to read file metadata: {}", input.display()))?;

    let file_size = metadata.len();
    let bytes_per_frame = num_sensors * 3 * bytes_per_sample;
    let total_samples = file_size as usize / bytes_per_frame;
    let total_length = total_samples as f64 / SAMPLE_RATE_HZ * PIG_SPEED_M_S;

    println!("  File Information");
    println!("  ──────────────────────────────────────────");
    println!("  Path:             {}", input.display());
    println!("  File size:        {:.2} MB", file_size as f64 / (1024.0 * 1024.0));
    println!("  Sensors:          {}", num_sensors);
    println!("  Bytes/sample:     {}", bytes_per_sample);
    println!("  Bytes/frame:      {}", bytes_per_frame);
    println!("  Total samples:    {}", total_samples);
    println!("  Sample rate:      {} Hz", SAMPLE_RATE_HZ);
    println!("  PIG speed:        {} m/s", PIG_SPEED_M_S);
    println!("  Pipeline length:  {:.2} m", total_length);
    println!("  Wall thickness:   {:.1} mm", WALL_THICKNESS_M * 1000.0);
    println!();

    Ok(())
}

fn print_banner(color_enabled: bool) {
    if color_enabled {
        println!();
        println!("  \x1b[38;5;39m╔══════════════════════════════════════════════════════╗\x1b[0m");
        println!("  \x1b[38;5;39m║\x1b[0m    \x1b[38;5;45m███╗   ███╗███████╗██╗         ██╗███╗   ██╗██╗   ██╗\x1b[0m   \x1b[38;5;39m║\x1b[0m");
        println!("  \x1b[38;5;39m║\x1b[0m    \x1b[38;5;45m████╗ ████║██╔════╝██║         ██║████╗  ██║██║   ██║\x1b[0m   \x1b[38;5;39m║\x1b[0m");
        println!("  \x1b[38;5;39m║\x1b[0m    \x1b[38;5;45m██╔████╔██║█████╗  ██║         ██║██╔██╗ ██║██║   ██║\x1b[0m   \x1b[38;5;39m║\x1b[0m");
        println!("  \x1b[38;5;39m║\x1b[0m    \x1b[38;5;45m██║╚██╔╝██║██╔══╝  ██║         ██║██║╚██╗██║██║   ██║\x1b[0m   \x1b[38;5;39m║\x1b[0m");
        println!("  \x1b[38;5;39m║\x1b[0m    \x1b[38;5;45m██║ ╚═╝ ██║███████╗███████╗    ██║██║ ╚████║╚██████╔╝\x1b[0m   \x1b[38;5;39m║\x1b[0m");
        println!("  \x1b[38;5;39m║\x1b[0m    \x1b[38;5;45m╚═╝     ╚═╝╚══════╝╚══════╝    ╚═╝╚═╝  ╚═══╝ ╚═════╝ \x1b[0m   \x1b[38;5;39m║\x1b[0m");
        println!("  \x1b[38;5;39m╚══════════════════════════════════════════════════════╝\x1b[0m");
        println!();
        println!("  \x1b[38;5;248mMFL Magnetic Dipole Inversion - Pipeline Defect Analysis\x1b[0m");
        println!("  \x1b[38;5;248mOffshore Oil & Gas Pipeline Inspection System\x1b[0m");
    } else {
        println!();
        println!("  ========================================================");
        println!("    MFL-INV  -  Magnetic Flux Leakage Inversion Tool");
        println!("    Pipeline Defect Analysis System");
        println!("  ========================================================");
    }
    println!();
}

fn atty_stdout() -> bool {
    #[cfg(target_family = "windows")]
    {
        false
    }
    #[cfg(not(target_family = "windows"))]
    {
        atty::is(atty::Stream::Stdout)
    }
}
