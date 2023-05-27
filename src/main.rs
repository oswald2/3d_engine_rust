extern crate sdl2;

use sdl2::event::Event;
use sdl2::gfx::primitives::DrawRenderer;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use sdl2::render::Canvas;
use sdl2::video::Window;
use std::f64::consts::PI;
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

    let camera = Vec3d::default();

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

            let line1 = Vec3d {
                x: tri_translated.p[1].x - tri_translated.p[0].x,
                y: tri_translated.p[1].y - tri_translated.p[0].y,
                z: tri_translated.p[1].z - tri_translated.p[0].z,
            };

            let line2 = Vec3d {
                x: tri_translated.p[2].x - tri_translated.p[0].x,
                y: tri_translated.p[2].y - tri_translated.p[0].y,
                z: tri_translated.p[2].z - tri_translated.p[0].z,
            };

            let mut normal = Vec3d {
                x: line1.y * line2.z - line1.z * line2.y,
                y: line1.z * line2.x - line1.x * line2.z,
                z: line1.x * line2.y - line1.y * line2.x,
            };

            let l = (normal.x * normal.x + normal.y * normal.y + normal.z * normal.z).sqrt();
            normal.x /= l;
            normal.y /= l;
            normal.z /= l;

            let cull = normal.x * (tri_translated.p[0].x - camera.x)
                + normal.y * (tri_translated.p[0].y - camera.y)
                + normal.z * (tri_translated.p[0].z - camera.z);

            if cull < 0.0 {
                let mut light_direction = Vec3d {
                    x: 0.0,
                    y: 0.0,
                    z: -1.0,
                };
                let l = (light_direction.x * light_direction.x
                    + light_direction.y * light_direction.y
                    + light_direction.z * light_direction.z)
                    .sqrt();
                light_direction.x /= l;
                light_direction.y /= l;
                light_direction.z /= l;

                let dp = normal.x * light_direction.x
                    + normal.y * light_direction.y
                    + normal.z * light_direction.z;

                let gs = ((dp + 1.0) * 0.5 * 255.0) as u8;
                //let gs = (dp * 255.0) as u8;
                let color = Color::RGB(gs, gs, gs);

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

                draw_triangle(&mut canvas, tri_projected, color);
            }
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

impl Default for Vec3d {
    fn default() -> Self {
        Vec3d {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
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
    canvas
        .filled_trigon(
            tr.p[0].x as i16,
            tr.p[0].y as i16,
            tr.p[1].x as i16,
            tr.p[1].y as i16,
            tr.p[2].x as i16,
            tr.p[2].y as i16,
            color,
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
