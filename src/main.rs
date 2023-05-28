extern crate sdl2;

use sdl2::event::Event;
use sdl2::gfx::primitives::DrawRenderer;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::render::Canvas;
use sdl2::video::Window;
use std::f64::consts::PI;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::num::{ParseFloatError, ParseIntError};
use std::ops::{Add, Mul, Sub};
use std::path::Path;
use std::time::Instant;

const WIDTH: u32 = 800;
const HEIGHT: u32 = 800;

pub fn main() {
    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();

    let window = video_subsystem
        .window("Rust 3D Engine", WIDTH, HEIGHT)
        .position_centered()
        .build()
        .unwrap();

    let mut canvas = window.into_canvas().build().unwrap();

    let mut event_pump = sdl_context.event_pump().unwrap();

    let mut camera = Vec3d::default();
    // let look_dir = Vec3d {
    //     x: 0.0,
    //     y: 0.0,
    //     z: 1.0,
    //     w: 1.0,
    // };
    let mut yaw = 0.0;

    let mut mesh = Mesh::new();
    match mesh.load_from_file(Path::new("axis.obj")) {
        Err(err) => {
            println!("Error loading file: {:?}", err);
            return;
        }
        _ => {}
    }

    let near = 0.1;
    let far = 1000.0;
    let aspect_ratio: f64 = HEIGHT as f64 / WIDTH as f64;

    let mat_proj = Mat4x4::make_projection(90.0, aspect_ratio, near, far);

    let mut theta: f64 = 0.0;

    let mut old = Instant::now();

    let mut vec_triangles_to_raster = Vec::new();

    'render: loop {
        let new = Instant::now();
        canvas.set_draw_color(Color::RGB(0, 0, 0));
        canvas.clear();

        let elapsed = (new - old).as_secs_f64();

        //theta += 1.0 * elapsed;

        let mat_rot_z = Mat4x4::make_rotation_z(theta * 0.5);
        let mat_rot_x = Mat4x4::make_rotation_x(theta);

        let mat_trans = Mat4x4::make_translation(0.0, 0.0, 8.0);

        let mut mat_world = &mat_rot_z * &mat_rot_x;
        mat_world = &mat_world * &mat_trans;

        let up = Vec3d {
            x: 0.0,
            y: 1.0,
            z: 0.0,
            w: 1.0,
        };
        let target = Vec3d {
            x: 0.0,
            y: 0.0,
            z: 1.0,
            w: 1.0,
        };
        let mat_camera_rot = Mat4x4::make_rotation_y(yaw);
        let look_dir = multiply_matrix_vector(&target, &mat_camera_rot);
        let target = camera + &look_dir;

        let mat_camera = Mat4x4::point_at(&camera, &target, &up);
        let mat_view = mat_camera.quick_inverse();

        for tri in &mesh.tris {
            let tri_transformed = Triangle {
                p: [
                    multiply_matrix_vector(&tri.p[0], &mat_world),
                    multiply_matrix_vector(&tri.p[1], &mat_world),
                    multiply_matrix_vector(&tri.p[2], &mat_world),
                ],
                color: Color::RGB(1, 1, 1),
            };

            let line1 = tri_transformed.p[1] - &tri_transformed.p[0];
            let line2 = tri_transformed.p[2] - &tri_transformed.p[0];

            let normal = line1.cross(&line2).normalise();

            let camera_ray = tri_transformed.p[0] - &camera;

            let cull = normal.dot(&camera_ray);

            if cull < 0.0 {
                let light_direction = Vec3d {
                    x: 0.0,
                    y: 0.0,
                    z: -1.0,
                    w: 1.0,
                }
                .normalise();

                let dp = normal.dot(&light_direction);
                let gs = ((dp + 1.0) * 0.5 * 255.0) as u8;
                let col = Color::RGB(gs, gs, gs);

                let tri_viewed = Triangle {
                    p: [
                        multiply_matrix_vector(&tri_transformed.p[0], &mat_view),
                        multiply_matrix_vector(&tri_transformed.p[1], &mat_view),
                        multiply_matrix_vector(&tri_transformed.p[2], &mat_view),
                    ],
                    color: col,
                };

                let mut tri_projected = Triangle {
                    p: [
                        multiply_matrix_vector(&tri_viewed.p[0], &mat_proj),
                        multiply_matrix_vector(&tri_viewed.p[1], &mat_proj),
                        multiply_matrix_vector(&tri_viewed.p[2], &mat_proj),
                    ],
                    color: col,
                };

                tri_projected.p[0] = tri_projected.p[0].vec_div(tri_projected.p[0].w);
                tri_projected.p[1] = tri_projected.p[1].vec_div(tri_projected.p[1].w);
                tri_projected.p[2] = tri_projected.p[2].vec_div(tri_projected.p[2].w);

                tri_projected.p[0].x *= -1.0;
                tri_projected.p[1].x *= -1.0;
                tri_projected.p[2].x *= -1.0;
                tri_projected.p[0].y *= -1.0;
                tri_projected.p[1].y *= -1.0;
                tri_projected.p[2].y *= -1.0;

                let offset_view = Vec3d {
                    x: 1.0,
                    y: 1.0,
                    z: 0.0,
                    w: 1.0,
                };

                // Scale into view.
                tri_projected.p[0] = tri_projected.p[0] + &offset_view;
                tri_projected.p[1] = tri_projected.p[1] + &offset_view;
                tri_projected.p[2] = tri_projected.p[2] + &offset_view;

                tri_projected.p[0].x *= 0.5 * WIDTH as f64;
                tri_projected.p[0].y *= 0.5 * HEIGHT as f64;
                tri_projected.p[1].x *= 0.5 * WIDTH as f64;
                tri_projected.p[1].y *= 0.5 * HEIGHT as f64;
                tri_projected.p[2].x *= 0.5 * WIDTH as f64;
                tri_projected.p[2].y *= 0.5 * HEIGHT as f64;

                vec_triangles_to_raster.push(tri_projected);

                // draw_triangle(&mut canvas, tri_projected, color);
            }
        }

        let forward = look_dir.vec_mul(160.0 * elapsed);

        for event in event_pump.poll_iter() {
            match event {
                Event::Quit { .. }
                | Event::KeyDown {
                    keycode: Some(Keycode::Escape),
                    ..
                } => break 'render,
                Event::KeyDown {
                    keycode: Some(Keycode::Up),
                    ..
                } => camera.y += 160.0 * elapsed,
                Event::KeyDown {
                    keycode: Some(Keycode::Down),
                    ..
                } => camera.y -= 160.0 * elapsed,
                Event::KeyDown {
                    keycode: Some(Keycode::Left),
                    ..
                } => camera.x -= 160.0 * elapsed,
                Event::KeyDown {
                    keycode: Some(Keycode::Right),
                    ..
                } => camera.x += 160.0 * elapsed,
                Event::KeyDown {
                    keycode: Some(Keycode::A),
                    ..
                } => yaw -= 80.0 * elapsed,
                Event::KeyDown {
                    keycode: Some(Keycode::D),
                    ..
                } => yaw += 80.0 * elapsed,
                Event::KeyDown {
                    keycode: Some(Keycode::W),
                    ..
                } => camera = camera + &forward,
                Event::KeyDown {
                    keycode: Some(Keycode::S),
                    ..
                } => camera = camera - &forward,
                _ => {}
            }
        }

        // Sort triangles
        vec_triangles_to_raster.sort_by(|t1, t2| {
            let z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0;
            let z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0;

            z2.total_cmp(&z1)
        });

        for tr in &vec_triangles_to_raster {
            draw_triangle(&mut canvas, tr);
        }
        vec_triangles_to_raster.clear();

        canvas.present();

        old = new;
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Vec3d {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
}

impl Default for Vec3d {
    fn default() -> Self {
        Vec3d {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            w: 1.0,
        }
    }
}

impl Vec3d {
    pub fn vec_mul(&self, k: f64) -> Vec3d {
        Vec3d {
            x: self.x * k,
            y: self.y * k,
            z: self.z * k,
            w: 1.0,
        }
    }

    pub fn vec_div(&self, k: f64) -> Vec3d {
        Vec3d {
            x: self.x / k,
            y: self.y / k,
            z: self.z / k,
            w: 1.0,
        }
    }

    pub fn dot(&self, rhs: &Vec3d) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    pub fn len(&self) -> f64 {
        self.dot(self).sqrt()
    }

    pub fn normalise(&self) -> Vec3d {
        let l = self.len();
        Vec3d {
            x: self.x / l,
            y: self.y / l,
            z: self.z / l,
            w: 1.0,
        }
    }

    pub fn cross(&self, rhs: &Vec3d) -> Vec3d {
        Vec3d {
            x: self.y * rhs.z - self.z * rhs.y,
            y: self.z * rhs.x - self.x * rhs.z,
            z: self.x * rhs.y - self.y * rhs.x,
            w: 1.0,
        }
    }
}

impl Add<&Vec3d> for Vec3d {
    type Output = Vec3d;

    fn add(self, rhs: &Vec3d) -> Self::Output {
        Vec3d {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
            w: self.w + rhs.w,
        }
    }
}

impl Sub<&Vec3d> for Vec3d {
    type Output = Vec3d;

    fn sub(self, rhs: &Vec3d) -> Self::Output {
        Vec3d {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
            w: self.w - rhs.w,
        }
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Triangle {
    pub p: [Vec3d; 3],
    pub color: Color,
}

pub struct Mesh {
    tris: Vec<Triangle>,
}

impl Mesh {
    pub fn new() -> Mesh {
        Mesh { tris: Vec::new() }
    }

    pub fn load_from_file(&mut self, file_path: &Path) -> Result<(), ObjParseError> {
        let file = File::open(file_path)?;
        let reader = BufReader::new(file);
        let mut vertices: Vec<Vec3d> = Vec::new();
        let mut triangles: Vec<Triangle> = Vec::new();

        for line in reader.lines() {
            if let Ok(line) = line {
                let fields: Vec<&str> = line.split_whitespace().collect();
                if fields.is_empty() {
                    continue;
                }

                match fields[0] {
                    "v" => {
                        if fields.len() < 4 {
                            return Err(ObjParseError::InvalidData(String::from(
                                "Invalid vertex line",
                            )));
                        }
                        let x = fields[1].parse::<f64>()?;
                        let y = fields[2].parse::<f64>()?;
                        let z = fields[3].parse::<f64>()?;
                        vertices.push(Vec3d { x, y, z, w: 1.0 });
                    }
                    "f" => {
                        if fields.len() < 4 {
                            return Err(ObjParseError::InvalidData(String::from(
                                "Invalid face line",
                            )));
                        }
                        let v1 = fields[1].parse::<usize>()?;
                        let v2 = fields[2].parse::<usize>()?;
                        let v3 = fields[3].parse::<usize>()?;
                        triangles.push(Triangle {
                            p: [vertices[v1 - 1], vertices[v2 - 1], vertices[v3 - 1]],
                            color: Color::RGB(255, 255, 255),
                        });
                    }
                    _ => continue,
                }
            }
        }

        self.tris = triangles;

        Ok(())
    }
}

#[derive(Debug)]
pub enum ObjParseError {
    IoError(std::io::Error),
    ParseError(ParseFloatError),
    ParseErrorInt(ParseIntError),
    InvalidData(String),
}

impl From<std::io::Error> for ObjParseError {
    fn from(error: std::io::Error) -> Self {
        ObjParseError::IoError(error)
    }
}

impl From<ParseFloatError> for ObjParseError {
    fn from(error: ParseFloatError) -> Self {
        ObjParseError::ParseError(error)
    }
}

impl From<ParseIntError> for ObjParseError {
    fn from(error: ParseIntError) -> Self {
        ObjParseError::ParseErrorInt(error)
    }
}

pub struct Mat4x4 {
    m: [[f64; 4]; 4],
}

impl Default for Mat4x4 {
    fn default() -> Self {
        Mat4x4 { m: [[0.0; 4]; 4] }
    }
}

impl Mat4x4 {
    pub fn identity() -> Mat4x4 {
        Mat4x4 {
            m: [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ],
        }
    }

    pub fn make_rotation_x(theta: f64) -> Mat4x4 {
        let mut matrix = Mat4x4::default();

        let cosf = theta.cos();
        let sinf = theta.sin();

        matrix.m[0][0] = 1.0;
        matrix.m[1][1] = cosf;
        matrix.m[1][2] = sinf;
        matrix.m[2][1] = -sinf;
        matrix.m[2][2] = cosf;
        matrix.m[3][3] = 1.0;

        matrix
    }

    pub fn make_rotation_y(theta: f64) -> Mat4x4 {
        let mut matrix = Mat4x4::default();

        let cosf = theta.cos();
        let sinf = theta.sin();

        matrix.m[0][0] = cosf;
        matrix.m[0][2] = sinf;
        matrix.m[2][0] = -sinf;
        matrix.m[1][1] = 1.0;
        matrix.m[2][2] = cosf;
        matrix.m[3][3] = 1.0;

        matrix
    }

    pub fn make_rotation_z(theta: f64) -> Mat4x4 {
        let mut matrix = Mat4x4::default();

        let cosf = theta.cos();
        let sinf = theta.sin();

        matrix.m[0][0] = cosf;
        matrix.m[0][1] = sinf;
        matrix.m[1][0] = -sinf;
        matrix.m[1][1] = cosf;
        matrix.m[2][2] = 1.0;
        matrix.m[3][3] = 1.0;

        matrix
    }

    pub fn make_translation(x: f64, y: f64, z: f64) -> Mat4x4 {
        let mut matrix = Mat4x4::default();

        matrix.m[0][0] = 1.0;
        matrix.m[1][1] = 1.0;
        matrix.m[2][2] = 1.0;
        matrix.m[3][3] = 1.0;
        matrix.m[3][0] = x;
        matrix.m[3][1] = y;
        matrix.m[3][2] = z;

        matrix
    }

    pub fn make_projection(fov_degress: f64, aspect_ratio: f64, near: f64, far: f64) -> Mat4x4 {
        let fov_rad = 1.0 / (fov_degress * 0.5 / 180.0 * PI).tan();
        let far_near = far - near;

        let mut matrix = Mat4x4::default();

        matrix.m[0][0] = aspect_ratio * fov_rad;
        matrix.m[1][1] = fov_rad;
        matrix.m[2][2] = far / far_near;
        matrix.m[3][2] = (-far * near) / far_near;
        matrix.m[2][3] = 1.0;
        matrix.m[3][3] = 0.0;

        matrix
    }

    pub fn point_at(pos: &Vec3d, target: &Vec3d, up: &Vec3d) -> Mat4x4 {
        let new_forward = (*target - pos).normalise();

        let a = new_forward.vec_mul(up.dot(&new_forward));
        let new_up = (*up - &a).normalise();

        let new_right = new_up.cross(&new_forward);

        let mut matrix = Mat4x4::default();

        matrix.m[0][0] = new_right.x;
        matrix.m[0][1] = new_right.y;
        matrix.m[0][2] = new_right.z;
        matrix.m[0][3] = 0.0;

        matrix.m[1][0] = new_up.x;
        matrix.m[1][1] = new_up.y;
        matrix.m[1][2] = new_up.z;
        matrix.m[1][3] = 0.0;

        matrix.m[2][0] = new_forward.x;
        matrix.m[2][1] = new_forward.y;
        matrix.m[2][2] = new_forward.z;
        matrix.m[2][3] = 0.0;

        matrix.m[3][0] = pos.x;
        matrix.m[3][1] = pos.y;
        matrix.m[3][2] = pos.z;
        matrix.m[3][3] = 1.0;

        matrix
    }

    pub fn quick_inverse(&self) -> Mat4x4 {
        let mut matrix = Mat4x4::default();

        matrix.m[0][0] = self.m[0][0];
        matrix.m[0][1] = self.m[1][0];
        matrix.m[0][2] = self.m[2][0];
        matrix.m[0][3] = 0.0;

        matrix.m[1][0] = self.m[0][1];
        matrix.m[1][1] = self.m[1][1];
        matrix.m[1][2] = self.m[2][1];
        matrix.m[1][3] = 0.0;

        matrix.m[2][0] = self.m[0][2];
        matrix.m[2][1] = self.m[1][2];
        matrix.m[2][2] = self.m[2][2];
        matrix.m[2][3] = 0.0;

        matrix.m[3][0] = -(self.m[3][0] * matrix.m[0][0]
            + self.m[3][1] * matrix.m[1][0]
            + self.m[3][2] * matrix.m[2][0]);
        matrix.m[3][1] = -(self.m[3][0] * matrix.m[0][1]
            + self.m[3][1] * matrix.m[1][1]
            + self.m[3][2] * matrix.m[2][1]);
        matrix.m[3][2] = -(self.m[3][0] * matrix.m[0][2]
            + self.m[3][1] * matrix.m[1][2]
            + self.m[3][2] * matrix.m[2][2]);
        matrix.m[3][3] = 1.0;

        matrix
    }
}

impl Mul<&Mat4x4> for &Mat4x4 {
    type Output = Mat4x4;

    fn mul(self, rhs: &Mat4x4) -> Self::Output {
        let mut matrix = Mat4x4::default();
        for c in 0..4 {
            for r in 0..4 {
                matrix.m[r][c] = self.m[r][0] * rhs.m[0][c]
                    + self.m[r][1] * rhs.m[1][c]
                    + self.m[r][2] * rhs.m[2][c]
                    + self.m[r][3] * rhs.m[3][c];
            }
        }
        matrix
    }
}

pub fn multiply_matrix_vector(i: &Vec3d, m: &Mat4x4) -> Vec3d {
    let x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
    let y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
    let z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
    let w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];

    Vec3d { x, y, z, w }
}

pub fn draw_triangle(canvas: &mut Canvas<Window>, tr: &Triangle) {
    canvas
        .filled_trigon(
            tr.p[0].x as i16,
            tr.p[0].y as i16,
            tr.p[1].x as i16,
            tr.p[1].y as i16,
            tr.p[2].x as i16,
            tr.p[2].y as i16,
            tr.color,
        )
        .unwrap();
    // canvas
    //     .trigon(
    //         tr.p[0].x as i16,
    //         tr.p[0].y as i16,
    //         tr.p[1].x as i16,
    //         tr.p[1].y as i16,
    //         tr.p[2].x as i16,
    //         tr.p[2].y as i16,
    //         Color::RGB(0, 0, 0),
    //     )
    //     .unwrap();
}
