extern crate sdl2;

use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::rect::Point;
use sdl2::render::Canvas;
use sdl2::video::Window;
use std::f64::consts::PI;
use std::time::{Duration, Instant};


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

    let mesh_cube = Mesh {
        tris: vec![
            // SOUTH
            Triangle {
                p: [
                    Vec3d {
                        x: 0.0,
                        y: 0.0,
                        z: 0.0,
                    },
                    Vec3d {
                        x: 0.0,
                        y: 1.0,
                        z: 0.0,
                    },
                    Vec3d {
                        x: 1.0,
                        y: 1.0,
                        z: 0.0,
                    },
                ],
            },
            Triangle {
                p: [
                    Vec3d {
                        x: 0.0,
                        y: 0.0,
                        z: 0.0,
                    },
                    Vec3d {
                        x: 1.0,
                        y: 1.0,
                        z: 0.0,
                    },
                    Vec3d {
                        x: 1.0,
                        y: 0.0,
                        z: 0.0,
                    },
                ],
            },
            // EAST
            Triangle {
                p: [
                    Vec3d {
                        x: 1.0,
                        y: 0.0,
                        z: 0.0,
                    },
                    Vec3d {
                        x: 1.0,
                        y: 1.0,
                        z: 0.0,
                    },
                    Vec3d {
                        x: 1.0,
                        y: 1.0,
                        z: 1.0,
                    },
                ],
            },
            Triangle {
                p: [
                    Vec3d {
                        x: 1.0,
                        y: 0.0,
                        z: 0.0,
                    },
                    Vec3d {
                        x: 1.0,
                        y: 1.0,
                        z: 1.0,
                    },
                    Vec3d {
                        x: 1.0,
                        y: 0.0,
                        z: 1.0,
                    },
                ],
            },
            // NORTH
            Triangle {
                p: [
                    Vec3d {
                        x: 1.0,
                        y: 0.0,
                        z: 1.0,
                    },
                    Vec3d {
                        x: 1.0,
                        y: 1.0,
                        z: 1.0,
                    },
                    Vec3d {
                        x: 0.0,
                        y: 1.0,
                        z: 1.0,
                    },
                ],
            },
            Triangle {
                p: [
                    Vec3d {
                        x: 1.0,
                        y: 0.0,
                        z: 1.0,
                    },
                    Vec3d {
                        x: 0.0,
                        y: 1.0,
                        z: 1.0,
                    },
                    Vec3d {
                        x: 0.0,
                        y: 0.0,
                        z: 1.0,
                    },
                ],
            },
            // WEST
            Triangle {
                p: [
                    Vec3d {
                        x: 0.0,
                        y: 0.0,
                        z: 1.0,
                    },
                    Vec3d {
                        x: 0.0,
                        y: 1.0,
                        z: 1.0,
                    },
                    Vec3d {
                        x: 0.0,
                        y: 1.0,
                        z: 0.0,
                    },
                ],
            },
            Triangle {
                p: [
                    Vec3d {
                        x: 0.0,
                        y: 0.0,
                        z: 1.0,
                    },
                    Vec3d {
                        x: 0.0,
                        y: 1.0,
                        z: 0.0,
                    },
                    Vec3d {
                        x: 0.0,
                        y: 0.0,
                        z: 0.0,
                    },
                ],
            },
            // TOP
            Triangle {
                p: [
                    Vec3d {
                        x: 0.0,
                        y: 1.0,
                        z: 0.0,
                    },
                    Vec3d {
                        x: 0.0,
                        y: 1.0,
                        z: 1.0,
                    },
                    Vec3d {
                        x: 1.0,
                        y: 1.0,
                        z: 1.0,
                    },
                ],
            },
            Triangle {
                p: [
                    Vec3d {
                        x: 0.0,
                        y: 1.0,
                        z: 0.0,
                    },
                    Vec3d {
                        x: 1.0,
                        y: 1.0,
                        z: 1.0,
                    },
                    Vec3d {
                        x: 1.0,
                        y: 1.0,
                        z: 0.0,
                    },
                ],
            },
            // BOTTOM
            Triangle {
                p: [
                    Vec3d {
                        x: 1.0,
                        y: 0.0,
                        z: 1.0,
                    },
                    Vec3d {
                        x: 0.0,
                        y: 0.0,
                        z: 1.0,
                    },
                    Vec3d {
                        x: 0.0,
                        y: 0.0,
                        z: 0.0,
                    },
                ],
            },
            Triangle {
                p: [
                    Vec3d {
                        x: 1.0,
                        y: 0.0,
                        z: 1.0,
                    },
                    Vec3d {
                        x: 0.0,
                        y: 0.0,
                        z: 0.0,
                    },
                    Vec3d {
                        x: 1.0,
                        y: 0.0,
                        z: 0.0,
                    },
                ],
            },
        ],
    };

    let near = 0.1;
    let far = 1000.0;
    let f_fov = 90.0;
    let aspect_ratio: f64 = HEIGHT as f64 / WIDTH as f64;
    let fov_rad: f64 = 1.0 / (f_fov * 0.5 / 180.0 * PI).tan() as f64;

    let mut mat_proj = Mat4x4 { m: [[0.0; 4]; 4] };

    mat_proj.m[0][0] = aspect_ratio * fov_rad;
    mat_proj.m[1][1] = fov_rad;
    mat_proj.m[2][2] = far / (far - near);
    mat_proj.m[3][2] = (-far * near) / (far - near);
    mat_proj.m[2][3] = 1.0;
    mat_proj.m[3][3] = 0.0;

    let mut theta: f64 = 0.0;

    let mut old = Instant::now();

    'render: loop {
        let new = Instant::now();

        canvas.set_draw_color(Color::RGB(0, 0, 0));
        canvas.clear();

        let elapsed = (new - old).as_secs_f64();

        theta += 1.0 * elapsed;

        let mut mat_rot_z = Mat4x4 { m: [[0.0; 4]; 4] };
        let mut mat_rot_x = Mat4x4 { m: [[0.0; 4]; 4] };

        mat_rot_z.m[0][0] = theta.cos();
        mat_rot_z.m[0][1] = theta.sin();
        mat_rot_z.m[1][0] = -(theta.sin());
        mat_rot_z.m[1][1] = theta.cos();
        mat_rot_z.m[2][2] = 1.0;
        mat_rot_z.m[3][3] = 1.0;

        mat_rot_x.m[0][0] = 1.0;
        mat_rot_x.m[1][1] = (theta * 0.5).cos();
        mat_rot_x.m[1][2] = (theta * 0.5).sin();
        mat_rot_x.m[2][1] = -((theta * 0.5).sin());
        mat_rot_x.m[2][2] = (theta * 0.5).cos();
        mat_rot_x.m[3][3] = 1.0;

        for tri in &mesh_cube.tris {
            let tri_rotated_z = Triangle {
                p: [
                    multiply_matrix_vector(&tri.p[0], &mat_rot_z),
                    multiply_matrix_vector(&tri.p[1], &mat_rot_z),
                    multiply_matrix_vector(&tri.p[2], &mat_rot_z),
                ],
            };

            let tri_rotated_zx = Triangle {
                p: [
                    multiply_matrix_vector(&tri_rotated_z.p[0], &mat_rot_x),
                    multiply_matrix_vector(&tri_rotated_z.p[1], &mat_rot_x),
                    multiply_matrix_vector(&tri_rotated_z.p[2], &mat_rot_x),
                ],
            };

            let mut tri_translated = tri_rotated_zx;

            tri_translated.p[0].z = tri_rotated_zx.p[0].z + 3.0;
            tri_translated.p[1].z = tri_rotated_zx.p[1].z + 3.0;
            tri_translated.p[2].z = tri_rotated_zx.p[2].z + 3.0;

            let mut tri_projected = Triangle {
                p: [
                    multiply_matrix_vector(&tri_translated.p[0], &mat_proj),
                    multiply_matrix_vector(&tri_translated.p[1], &mat_proj),
                    multiply_matrix_vector(&tri_translated.p[2], &mat_proj),
                ],
            };

            // Scale into view.
            tri_projected.p[0].x += 1.0;
            tri_projected.p[0].y += 1.0;
            tri_projected.p[1].x += 1.0;
            tri_projected.p[1].y += 1.0;
            tri_projected.p[2].x += 1.0;
            tri_projected.p[2].y += 1.0;

            tri_projected.p[0].x *= 0.5 * WIDTH as f64;
            tri_projected.p[0].y *= 0.5 * HEIGHT as f64;
            tri_projected.p[1].x *= 0.5 * WIDTH as f64;
            tri_projected.p[1].y *= 0.5 * HEIGHT as f64;
            tri_projected.p[2].x *= 0.5 * WIDTH as f64;
            tri_projected.p[2].y *= 0.5 * HEIGHT as f64;

            draw_triangle(&mut canvas, tri_projected, Color::RGB(255, 255, 255));
        }

        for event in event_pump.poll_iter() {
            match event {
                Event::Quit { .. }
                | Event::KeyDown {
                    keycode: Some(Keycode::Escape),
                    ..
                } => break 'render,
                //Event::KeyDown { keycode: Some(Keycode::A), .. } => theta += 0.1,
                _ => {}
            }
        }
        canvas.present();

        old = new;
    }
}

#[derive(Debug, Copy, Clone)]
pub struct Vec3d {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Debug, Copy, Clone)]
pub struct Triangle {
    pub p: [Vec3d; 3],
}

pub struct Mesh {
    tris: Vec<Triangle>,
}

pub struct Mat4x4 {
    m: [[f64; 4]; 4],
}

pub fn multiply_matrix_vector(i: &Vec3d, m: &Mat4x4) -> Vec3d {
    let x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
    let y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
    let z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];
    let w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];

    if w != 0.0 {
        Vec3d {
            x: x / w,
            y: y / w,
            z: z / w,
        }
    } else {
        Vec3d { x, y, z }
    }
}

pub fn draw_triangle(canvas: &mut Canvas<Window>, tr: Triangle, color: Color) {
    let point1 = Point::new(tr.p[0].x as i32, tr.p[0].y as i32);
    let point2 = Point::new(tr.p[1].x as i32, tr.p[1].y as i32);
    let point3 = Point::new(tr.p[2].x as i32, tr.p[2].y as i32);

    canvas.set_draw_color(color);

    canvas.draw_line(point1, point2).unwrap();
    canvas.draw_line(point2, point3).unwrap();
    canvas.draw_line(point3, point1).unwrap();
}
