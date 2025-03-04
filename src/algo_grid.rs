use alloc::collections::VecDeque;
use alloc::string::{String, ToString};
use alloc::vec::Vec;
use core::iter::repeat;
use core::ops::SubAssign;

use rayon::prelude::*;

use base16ct;
use cgmath::MetricSpace;
#[allow(unused_imports)]
use cgmath::num_traits::Float;
use cgmath::Point2;
use ndarray::{Array1, Array2, Array3, ArrayBase, ArrayView1, ArrayView2, Axis};
use ndarray::arr1;
use peroxide::fuga::*;
use sha2::{Digest, Sha256};
use tri_mesh::mesh::Mesh;

use crate::contour::Rect;
use crate::polyline::GenPolyLines;

type VectorTriangles = Array3<f64>;
//type Triangle = Array3<f64>;
// [[x0, y0, z0], [x1, y1, z1], [x3, y3, z3]]
// [[x0, y0, z0], [x1, y1, z1], [x3, y3, z3]]
// [[x0, y0, z0], [x1, y1, z1], [x3, y3, z3]]
type Vec2 = Point2<f64>;


pub(crate) fn find_top_std(centers: &Vec<Vec<Vec2>>, depth: usize, grid_size: i16, rect: Rect) -> Vec<String> {
    let mut hashes = vec![];
    if centers.len() == 0 {
        return hashes;
    }

    let ss = GenPolyLines::select_top(centers, depth, grid_size, rect);

    for a in ss.iter() {
        let data: Vec<u8> = a.1.nodes.as_slice().iter()
            .flat_map(|&p| [p.x.to_be_bytes(), p.y.to_be_bytes()])
            .flatten()
            .collect();

        let mut hasher = Sha256::new();
        hasher.update(data.as_slice());

        let mut buf = [0u8; 64];
        let hash = hasher.finalize();
        let hex_hash = base16ct::lower::encode_str(&hash, &mut buf).unwrap();

        hashes.push(hex_hash.to_string());
    }
    hashes.dedup();
    hashes
}

pub(crate) fn find_top_std_2(centers: &Vec<Vec<Vec2>>, depth: usize, n_sect: usize, grid_size: usize, rect: Rect) -> Vec<String> {
    let mut hashes = vec![];
    if centers.len() == 0 {
        return hashes;
    }

    const N: usize = 2;
    let ss = GenPolyLines::select_top_all(centers, depth, grid_size, rect);
    if ss.len() < n_sect {
        return hashes;
    }
    let mut best_totals: VecDeque<(f64, Vec<u8>)> = VecDeque::with_capacity(depth);

    let mut ff = |d: f64, hash: Vec<u8>| {
        let len = best_totals.len();
        if len > 0 {
            if d < best_totals.get(len - 1).unwrap().0 || len <= depth {
                if len == depth {
                    best_totals.pop_front();
                }
                best_totals.push_back((d, hash));
                best_totals.make_contiguous()
                    .sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
            }
        } else {
            best_totals.push_back((d, hash));
        }
    };

    let mut stack: Vec<usize> = repeat(0).take(n_sect).collect();

    loop {
        let mut sco = 0.;
        let mut h: Vec<u8> = Vec::new();
        for l in 0..n_sect {
            let k = stack[l];
            if k < ss[l].len() {
                sco += ss[l][k].0;
                h.extend(ss[l][k].1.clone());
            }
        }
        ff(sco, h);

        let mut j = 0;
        while j < n_sect {
            if stack[j] < N - 1 {
                stack[j] += 1;
                break;
            }
            stack[j] = 0;
            j += 1;
        }
        if j == n_sect {
            break;
        }
    }

    for hash in best_totals.iter() {
        let mut hasher = Sha256::new();
        hasher.update(hash.1.as_slice());

        let mut a: Vec<u8> = repeat(0).take(32 * n_sect).collect();
        let mut buf = a.as_mut();
        let hash = hasher.finalize();
        let hex_hash = base16ct::lower::encode_str(&hash, &mut buf).unwrap();

        hashes.push(hex_hash.to_string());
    }
    hashes.dedup();
    hashes
}

pub(crate) fn find_top_std_3(centers: &Vec<Vec<Vec2>>, depth: usize, n_sect: usize, grid_size: usize, rect: Rect) -> Vec<String> {
    let mut hashes = vec![];
    if centers.len() == 0 {
        return hashes;
    }

    const N: usize = 2;
    let ss = GenPolyLines::select_top_all_3(centers, depth, grid_size, rect);
    if ss.len() < n_sect {
        return hashes;
    }

    let mut best_totals: Vec<(f64, Vec<u8>)> = Vec::with_capacity(depth);

    let mut ff = |d: f64, hash: Vec<u8>| {
        if let Some(_) = best_totals.iter().find(|a| a.0 == d) {
            return;
        } else {
            if best_totals.len() == depth {
                let m = best_totals.iter()
                    .enumerate()
                    .max_by(|(_, a), (_, b)|
                        a.0.partial_cmp(&b.0).unwrap_or(core::cmp::Ordering::Equal)
                    )
                    .map(|(index, _)| index);

                if let Some(i) = m {
                    best_totals[i] = (d, hash);
                }
            } else {
                best_totals.push((d, hash));
            }
        }
    };

    let mut stack: Vec<usize> = repeat(0).take(n_sect).collect();

    loop {
        let mut sco = 0.;
        let mut h: Vec<u8> = Vec::new();
        for l in 0..n_sect {
            let k = stack[l];
            if k < ss[l].len() {
                sco += ss[l][k].0;
                h.extend(ss[l][k].1.clone());
            }
        }
        ff(sco, h);

        let mut j = 0;
        while j < n_sect {
            if stack[j] < N - 1 {
                stack[j] += 1;
                break;
            }
            stack[j] = 0;
            j += 1;
        }
        if j == n_sect {
            break;
        }
    }

    best_totals.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    for hash in best_totals.iter() {
        let mut hasher = Sha256::new();
        hasher.update(hash.1.as_slice());

        let mut a: Vec<u8> = repeat(0).take(32 * n_sect).collect();
        let mut buf = a.as_mut();
        let hash = hasher.finalize();
        let hex_hash = base16ct::lower::encode_str(&hash, &mut buf).unwrap();

        hashes.push(hex_hash.to_string());
    }
    hashes.dedup();
    hashes
}

pub(crate) fn find_top_std_4(
    cntrs: &Vec<Vec<Vec2>>, depth: usize, n_sect: usize, grid_size: usize, rect: Rect,
) -> Vec<String> {
    let mut hashes = vec![];
    if cntrs.len() == 0 {
        return hashes;
    }

    const N: usize = 2;
    let ss = GenPolyLines::select_top_all_4(cntrs, N, grid_size, rect);
    if ss.len() < n_sect {
        return hashes;
    }

    let mut best_totals: Vec<(f64, Vec<u8>)> = Vec::with_capacity(depth);

    let mut ff = |d: f64, hash: Vec<u8>| {
        if let Some(_) = best_totals.iter().find(|a| a.0 == d) {
            return
        }
        else {
            if best_totals.len() == depth {
                let m = best_totals.iter()
                    .enumerate()
                    .max_by(|(_, a), (_, b)|
                        a.0.partial_cmp(&b.0).unwrap_or(core::cmp::Ordering::Equal)
                    );

                if let Some((i, r)) = m {
                    if r.0 > d {
                        best_totals[i] = (d, hash);
                    }
                }
            } else {
                best_totals.push((d, hash));
            }
        }
    };

    let mut stack: Vec<usize> = repeat(0).take(n_sect).collect();

    loop {
        let mut sco = 0.;
        let mut h: Vec<u8> = Vec::new();
        for l in 0..n_sect {
            let k = stack[l];
            if k < ss[l].len() {
                sco += ss[l][k].0;
                h.extend(ss[l][k].1.clone());
            }
        }
        ff(sco, h);

        let mut j = 0;
        while j < n_sect {
            if ss[j].len() == 0 {
                j += 1;
                continue
            }
            if stack[j] < N - 1 {
                stack[j] += 1;
                break
            }
            stack[j] = 0;
            j += 1;
        }
        if j == n_sect {
            break
        }
    }

    best_totals.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    for hash in best_totals.iter() {
        let mut hasher = Sha256::new();
        hasher.update(hash.1.as_slice());

        let mut a: Vec<u8> = repeat(0).take(32 * n_sect).collect();
        let mut buf= a.as_mut();
        let hash = hasher.finalize();
        let hex_hash = base16ct::lower::encode_str(&hash, &mut buf).unwrap();

        hashes.push(hex_hash.to_string());
    }
    hashes.dedup();
    hashes
}


fn cross(triangles: &VectorTriangles) -> Array2<f64> {
    let dims = triangles.dim();
    let mut d = Array3::zeros((dims.0, 2, dims.1));

    for (i, m) in triangles.axis_iter(Axis(0)).enumerate() {
        for (j, p) in m.axis_iter(Axis(1)).enumerate() {
            d[[i, 0, j]] = p[1] - p[0];
            d[[i, 1, j]] = p[2] - p[1];
        }
    }

    let mut cr = Array2::zeros((dims.0, dims.1));
    for (i, m) in d.axis_iter(Axis(0)).enumerate() {
        let a: ArrayView1<f64> = m.slice(s![0, ..]);
        let b: ArrayView1<f64> = m.slice(s![1, ..]);
        cr[[i, 0]] = a[1] * b[2] - a[2] * b[1];
        cr[[i, 1]] = a[2] * b[0] - a[0] * b[2];
        cr[[i, 2]] = a[0] * b[1] - a[1] * b[0];
    }
    cr
}

pub fn mass_properties(triangles: VectorTriangles) -> (Array1<f64>, Array2<f64>) {

    let p0: ArrayView2<f64> = triangles.slice(s![0.., 0, 0..]);
    let p1: ArrayView2<f64> = triangles.slice(s![0.., 1, 0..]);
    let p2: ArrayView2<f64> = triangles.slice(s![0.., 2, 0..]);

    let f1 = &p0 + &p1 + &p2;
    let f2 = &p0 * &p0 + &p1 * &p1 + &p0 * &p1 + &p1 * &f1;
    let f3 = &p0 * &p0 * &p0 + &p0 * &p0 * &p1 + &p0 * &p1 * &p1 + &p1 * &p1 * &p1 + &p2 * &f2;

    let g0 = &(&f2 + &(&(&p0 + &f1) * &p0));
    let g1 = &(&f2 + &(&(&p1 + &f1) * &p1));
    let g2 = &(&f2 + &(&(&p2 + &f1) * &p2));

    let d = f1.nrows();
    let mut integral: Array2<f64> = Array2::zeros((10, d));

    let crosses = cross(&triangles);

    integral.slice_mut(s![0..1, ..]).assign(&(&crosses.slice(s![.., 0]) * &f1.slice(s![.., 0])));
    integral.slice_mut(s![1..4, ..]).assign(&(&crosses * &f2).t().slice(s![.., ..]));
    integral.slice_mut(s![4..7, ..]).assign(&(&crosses * &f3).t().slice(s![.., ..]));

    let results: Vec<Array2<f64>> = (0..3).into_par_iter().map(|i| {
        let triangle_i = (i + 1) % 3;
        let mut local_integral: Array2<f64> = Array2::zeros((10, d));
        local_integral.slice_mut(s![i+7, ..]).assign(
            &(&crosses.slice(s![.., i]) * &(
                &triangles.slice(s![.., 0, triangle_i]) * &g0.slice(s![.., i]) +
                &triangles.slice(s![.., 0, triangle_i]) * &g1.slice(s![.., i]) +
                &triangles.slice(s![.., 0, triangle_i]) * &g2.slice(s![.., i]))
            )
        );
        local_integral
    }).collect();

    for result in results {
        integral += &result;
    }

    let coefficients: Array1<f64> = arr1(&[1. / 6., 1. / 24., 1. / 24., 1. / 24., 1. / 60., 1. / 60., 1. / 60., 1. / 120., 1. / 120., 1. / 120.]);
    let integrated: Array1<f64> = integral.sum_axis(Axis(1)) * coefficients;
    let volume = integrated[0];
    let center_mass: Array1<f64> = if volume.abs() < 1e-10 {
        arr1(&[0., 0., 0.])
    } else {
        let a = &integrated.slice(s![1..4]);
        a / volume
    };

    let density = 1.0;

    let mut inertia: Array2<f64> = Array2::zeros((3, 3));

    inertia[[0, 0]] = integrated[5] + integrated[6] - volume * (center_mass[1].powi(2) + center_mass[2].powi(2));
    inertia[[1, 1]] = integrated[4] + integrated[6] - volume * (center_mass[0].powi(2) + center_mass[2].powi(2));
    inertia[[2, 2]] = integrated[4] + integrated[5] - volume * (center_mass[0].powi(2) + center_mass[1].powi(2));
    inertia[[0, 1]] = integrated[7] - volume * center_mass[0] * center_mass[1];
    inertia[[1, 2]] = integrated[8] - volume * center_mass[1] * center_mass[2];
    inertia[[0, 2]] = integrated[9] - volume * center_mass[0] * center_mass[2];
    inertia[[2, 0]] = inertia[[0, 2]];
    inertia[[2, 1]] = inertia[[1, 2]];
    inertia[[1, 0]] = inertia[[0, 1]];
    inertia *= density;

    (center_mass, inertia)
}



fn principal_axis(inertia: Array2<f64>) -> (Array1<f64>, Array2<f64>) {
    let negate_non_diagonal: Array2<f64> = &(Array2::eye(3) * 2.0) - 1.0;
    let a: Array2<f64> = inertia * negate_non_diagonal;

    let m = matrix(a.as_slice().unwrap().to_vec(), 3, 3, Row);
    let e = eigen(&m, Jacobi);
    let (c, v) = e.extract();

    let components = arr1(c.as_slice());
    let vectors: Array2<f64> = ArrayBase::from_shape_vec((3, 3), v.data).unwrap();

    // eigen returns them as column vectors, change them to row vectors
    (components, vectors.reversed_axes())
}

#[allow(dead_code)]
fn transform_around(matrix: Array2<f64>, point: &Array1<f64>) -> Array2<f64> {
    let mut translate: Array2<f64> = Array2::eye(4);
    translate.slice_mut(s![..3, ..3]).sub_assign(point);

    let mut result = matrix.dot(&translate);
    translate.slice_mut(s![..3, 3]).assign(point);
    result = translate.dot(&result);

    return result;
}

pub fn principal_inertia_transform(triangles: VectorTriangles) -> Array2<f64> {
    let (center_mass, inertia) = mass_properties(triangles);
    let (_components, vectors) = principal_axis(inertia);

    // TODO: Reorder vectors by components
    let mut transform = Array2::eye(4);

    // TODO:
    transform.slice_mut(s![..3, ..3]).assign(&vectors);

    // let mut tr = transform_around(transform, &center_mass);
    transform.slice_mut(s![..3, 3]).sub_assign(&center_mass);

    transform
}

pub fn intersect(mesh: &Mesh, z_sect: f64) -> Vec::<Vec2> {
    let mut sect = Vec::<Vec2>::new();

    for vertex_id in mesh.vertex_iter() {
        let p = mesh.vertex_position(vertex_id);
        if (p.z - z_sect).abs() < 0.15 {
            sect.push(Vec2{x: p.x, y: p.y});
        }
    }
    sect
}

pub fn intersect_2(mesh: &Mesh, z_sect: f64, delta: f64) -> Vec::<Vec2> {
    let mut sect = Vec::<Vec2>::new();

    for edge_id in mesh.edge_iter() {
        let (p1, p2) = mesh.edge_positions(edge_id);
        if p2.z >= z_sect && p1.z <= z_sect || p2.z <= z_sect && p1.z >= z_sect {
            let (x, y);
            let z1 = z_sect - p1.z;
            let z2 = p2.z - z_sect;
            if z1.abs() < delta {
                (x, y) = (p1.x, p1.y);
            }
            else if z2.abs() < delta {
                (x, y) = (p2.x, p2.y);
            }
            else {
                let k = z2 / z1;
                x = (p2.x + k * p1.x) / (k + 1.0);
                y = (p2.y + k * p1.y) / (k + 1.0);
            }
            sect.push(Vec2{x, y});
        }
    }
    sect
}


pub fn get_contour(sect: Vec<Vec2>) -> Vec<Point2<f64>> {
    let len = sect.len();

    // Compute mt in parallel
    let mt: Vec<Vec<f64>> = sect.par_iter().enumerate().map(|(i, p)| {
        let mut local_mt = vec![0f64; len];
        for (j, q) in sect.iter().enumerate() {
            local_mt[j] = p.distance2(*q);
        }
        local_mt
    }).collect();

    let mut ii: Vec<usize> = (0..len).collect();
    ii.par_sort_unstable_by(|&i, &j| {
        let v_i = &mt[i];
        let v_j = &mt[j];
        v_i[j].partial_cmp(&v_j[i]).unwrap()
    });

    let mut cntr: Vec<Point2<f64>> = sect
        .par_iter().enumerate()
        .map(|(i, &_a)| sect[ii[i]])
        .collect();

    cntr.push(*cntr.first().unwrap());

    let p0 = *cntr.first().unwrap();
    let pn = *cntr.last().unwrap();
    let d = cntr.first().unwrap().distance(*cntr.last().unwrap()).sqrt();
    let d2 = cntr[0].distance(cntr[1]).sqrt();

    let nn = (d / d2) as i32;
    cntr.par_extend((0..nn).into_par_iter().map(|n| {
        let k = (pn.y - p0.y) / (pn.x - p0.x);
        Point2 { x: p0.x + (n as f64) * d2, y: p0.y + (n as f64) * d2 * k }
    }));

    cntr
}


