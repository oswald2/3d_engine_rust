extern crate sdl2;

use sdl2::event::Event;
use sdl2::keyboard::Keycode;
use sdl2::pixels::Color;
use std::i128::MAX;
use std::time::Duration;

use num::complex::Complex;
use rand::prelude::*;


const WIDTH: u32 = 1280;
const HEIGHT: u32 = 1024;
const MAX_ITER: u32 = 1024;

pub fn main() {
    let sdl_context = sdl2::init().unwrap();
    let video_subsystem = sdl_context.video().unwrap();

    let window = video_subsystem
        .window("rust-sdl2 demo", WIDTH, HEIGHT)
        .position_centered()
        .build()
        .unwrap();

    let mut canvas = window.into_canvas().build().unwrap();

    canvas.set_draw_color(Color::RGB(0, 0, 0));
    canvas.clear();
    canvas.present();

    let mut event_pump = sdl_context.event_pump().unwrap();

    let mut rng = rand::thread_rng();
    let mut colors = Vec::new();
    for i in 0..MAX_ITER {
        let r: u8 = rng.gen::<u8>();
        let g: u8 = rng.gen::<u8>();
        let b: u8 = rng.gen::<u8>();
        colors.push(Color::RGB(r, g, b));
    }

    for px in 0..WIDTH {
        for py in 0..HEIGHT {
            let x0 = pix_to_coord_x(px, 0, WIDTH, -2.0, 1.0);
            let y0 = pix_to_cord_y(py, 0, HEIGHT, -1.5, 1.5);

            let mut x = 0.0;
            let mut y = 0.0;

            let mut iter = 0;

            let mut x2 = 0.0;
            let mut y2 = 0.0;

            while (x2 + y2 <= 4.0) && (iter < MAX_ITER) {
                y = 2.0 * x * y + y0;
                x = x2 - y2 + x0;
                x2 = x * x;
                y2 = y * y;
                iter += 1;
            }

            if (iter >= MAX_ITER) {
                canvas.set_draw_color(Color::RGB(0, 0, 0));
            } else {
                // let val: u8 = (iter as f64 / MAX_ITER as f64 * 255.0) as u8;
                // canvas.set_draw_color(Color::RGB(val, val, val));
                canvas.set_draw_color(colors[iter as usize]);
            }

            canvas.draw_point((px as i32, py as i32)).unwrap();
        }
    }

    'running: loop {
        for event in event_pump.poll_iter() {
            match event {
                Event::Quit { .. }
                | Event::KeyDown {
                    keycode: Some(Keycode::Escape),
                    ..
                } => break 'running,
                _ => {}
            }
        }
        canvas.present();
        ::std::thread::sleep(Duration::new(0, 1_000_000_000u32 / 60));
    }
}

fn pix_to_coord_x(x: u32, min: u32, max: u32, range_min: f64, range_max: f64) -> f64 {
    let rn = range_max - range_min;
    let ro = (max - min) as f64;

    (x - min) as f64 / ro * rn + range_min
}

fn pix_to_cord_y(y: u32, min: u32, max: u32, range_min: f64, range_max: f64) -> f64 {
    let rn = range_max - range_min;
    let ro = (max - min) as f64;

    (max - y) as f64 / ro * rn + range_min
}
