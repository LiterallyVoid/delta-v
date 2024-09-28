#![feature(try_blocks)]
#![feature(let_chains)]

type Semimajors<'a> = [(&'a str, f64, Option<&'a Planet<'a>>)];

type Hops<'a> = Vec<(&'a str, Vec<f64>)>;

trait Thousands {}

/// Gravitational constant, in units of N·m²·kg⁻².
const G: f64 = 6.674_30e-11;

// v² = µ(2/r - 1/a)

#[derive(Copy, Clone, Debug)]
pub struct Planet<'a> {
    pub name: &'a str,

    /// Orbit body and semi-major.
    pub orbit: Option<(&'a Planet<'a>, f64)>,

    /// This planet's mass in kilograms.
    pub mass: f64,

    /// The distance between this planet's center and its sea level on the equator, in meters.
    pub radius: f64,

    pub sidereal_day: f64,

    /// Low orbit altitude, in meters.
    pub low_orbit: f64,
}

impl Planet<'_> {
    fn low_orbit_velocity(&self) -> f64 {
        let r = self.radius + self.low_orbit;
        let mu = G * self.mass;

        (mu * (2.0 / r - 1.0 / r)).sqrt()
    }

    fn escape_velocity(&self) -> f64 {
        let r = self.radius + self.low_orbit;
        let mu = G * self.mass;

        (mu * (2.0 / r)).sqrt()
    }
}

#[derive(Copy, Clone, Debug)]
pub struct Orbit<'a> {
    pub planet: &'a Planet<'a>,

    // The periapsis of this orbit, in meters.
    pub periapsis: f64,

    // The apoapsis of this orbit, in meters.
    // This may be `infinity` for a parabolic escape trajectory.
    pub apoapsis: f64,
}

impl Orbit<'_> {
    pub fn semi_major(&self) -> f64 {
        (self.periapsis + self.apoapsis) / 2.0 + self.planet.radius
    }

    pub fn velocity_at_altitude(&self, altitude: f64) -> f64 {
        let r = altitude + self.planet.radius;
        (G * self.planet.mass * (2.0 / r - 1.0 / self.semi_major())).sqrt()
    }

    fn burn(&self, altitude: f64, new_orbit: Self) -> (f64, Self) {
        assert!(altitude >= self.periapsis);
        assert!(altitude <= self.apoapsis);

        (
            orbit_change_delta_v(
                self.planet,
                altitude,
                self.semi_major(),
                new_orbit.semi_major(),
            ),
            new_orbit,
        )
    }

    pub fn with_apoapsis(&self, apoapsis: f64) -> (f64, Self) {
        assert!(apoapsis >= self.periapsis);

        self.burn(self.periapsis, Orbit { apoapsis, ..*self })
    }

    /// BUG: this assumes `altitude` is higher than `apoapsis`.
    pub fn with_hohmann(&self, altitude: f64) -> (f64, Self) {
        let dv1;
        let dv2;
        let mut orbit;

        (dv1, orbit) = self.with_apoapsis(altitude);
        (dv2, orbit) = orbit.with_circularization();

        (dv1 + dv2, orbit)
    }

    pub fn with_hohmann_2(&self, altitude: f64) -> (f64, f64, Self) {
        let dv1;
        let dv2;
        let mut orbit;

        (dv1, orbit) = self.with_apoapsis(altitude);
        (dv2, orbit) = orbit.with_circularization();

        (dv1, dv2, orbit)
    }

    pub fn with_escape(&self) -> (f64, Self) {
        self.burn(
            self.periapsis,
            Orbit {
                apoapsis: f64::INFINITY,
                ..*self
            },
        )
    }

    /// Raise periapsis to apoapsis.
    pub fn with_circularization(&self) -> (f64, Self) {
        self.burn(
            self.apoapsis,
            Orbit {
                periapsis: self.apoapsis,
                ..*self
            },
        )
    }
}

/// Returns the ΔV required for an orbit change from a semi-major of `semi_major_initial` m to a
/// semi-major of `semi_major_final` m, where the impulse is executed at `altitude` m above the
/// planet's surface.
///
/// Semi-major can be calculated as the average of apoapsis and periapsis altitudes, plus the
/// planet's radius.
pub fn orbit_change_delta_v(
    planet: &Planet,
    altitude: f64,
    semi_major_initial: f64,
    semi_major_final: f64,
) -> f64 {
    let mu = G * planet.mass;
    let r = altitude + planet.radius;

    let x = 2.0 / r;

    let v_0 = (x - 1.0 / semi_major_initial).sqrt();
    let v_1 = (x - 1.0 / semi_major_final).sqrt();

    // sqrt(a + b) = x + y = sqrt(x^2 + y^2)

    mu.sqrt() * (v_1 - v_0).abs()
}

/// Return the nominal orbital velocity of a satellite that's orbiting around `planet`, with a
/// semi major of `semi_major`.
pub fn velocity_from_semi_major(planet: &Planet, semi_major: f64) -> f64 {
    let mu = G * planet.mass;
    (mu * (2.0 / semi_major - 1.0 / semi_major)).sqrt()
}

pub fn calculate_escape_dv(planet: &Planet, dv: f64) -> f64 {
    let mu = G * planet.mass;
    let r = planet.radius + planet.low_orbit;
    let periapsis_velocity = (dv.powf(2.0) + 2.0 * (mu / r)).sqrt();

    assert!(periapsis_velocity > planet.escape_velocity());

    periapsis_velocity - planet.escape_velocity()
}

fn main() {
    let sun = Planet {
        name: "Sun",
        orbit: None,

        mass: 1.988_5e+30,
        radius: 6.957e+8,

        sidereal_day: 86_400.0 * 25.0,
        low_orbit: 10_000_000.0,
    };

    let mercury = Planet {
        name: "Mercury",
        orbit: Some((&sun, 57_910_000_000.0)),
        mass: 3.3011e+23,
        radius: 2_439_700.0,
        sidereal_day: 86_400.0 * 58.65,
        low_orbit: 100_000.0,
    };

    let venus = Planet {
        name: "Venus",
        orbit: Some((&sun, 108_200_000_000.0)),
        mass: 4.8675e+24,
        radius: 6_051_800.0,
        sidereal_day: 86_400.0 * 243.0226,
        low_orbit: 250_000.0,
    };

    let earth = Planet {
        name: "Earth",
        orbit: Some((&sun, 149_598_023_000.0)),
        mass: 5.972_168e+24,
        radius: 6_378_137.0,
        sidereal_day: 86_164.090_500,
        low_orbit: 250_000.0,
    };

    let mars = Planet {
        name: "Mars",
        orbit: Some((&sun, 227_940_000_000.0)),
        mass: 6.4171e+23,
        radius: 3_389_500.0,
        sidereal_day: 86_400.0 * 1.025_927,
        low_orbit: 250_000.0,
    };

    let jupiter = Planet {
        name: "Jupiter",
        orbit: Some((&sun, 778_500_000_000.0)),
        mass: 1.8982e+27,
        radius: 69_911_000.0,
        sidereal_day: 3_600.0 * 9.9250,
        low_orbit: 8_000_000.0,
    };

    let saturn = Planet {
        name: "Saturn",
        orbit: Some((&sun, 1_429_000_000_000.0)),
        mass: 5.6834e+26,
        radius: 58_232_000.0,
        // not sidereal but a saturn year is like 30 earth years so
        // the effect will be negligible.
        sidereal_day: 3_600.0 * 10.7,
        low_orbit: 1_000_000.0,
    };

    let uranus = Planet {
        name: "Uranus",
        orbit: Some((&sun, 2_870_000_000_000.0)),
        mass: 8.6810e+25,
        radius: 25_362_000.0,
        sidereal_day: 3_600.0 * 17.0,
        low_orbit: 1_000_000.0,
    };

    let neptune = Planet {
        name: "Neptune",
        orbit: Some((&sun, 4_498_000_000_000.0)),
        mass: 1.02413e+26,
        radius: 24_622_000.0,
        sidereal_day: 3_600.0 * 16.0,
        low_orbit: 1_000_000.0,
    };

    let semimajors: &Semimajors = &[
        &mercury, &venus, &earth, &mars, &jupiter, &saturn, &uranus, &neptune,
    ]
    .into_iter()
    .map(|planet| (planet.name, planet.orbit.unwrap().1, Some(planet)))
    .collect::<Vec<_>>();

    fn launch<'a>(planet: &'a Planet<'a>) -> (f64, Orbit<'a>) {
        let orbit = Orbit {
            planet,
            periapsis: 0.0,
            apoapsis: 0.0,
        };
        let dv0 = orbit.velocity_at_altitude(orbit.periapsis);
        let (dv1, orbit) = orbit.with_apoapsis(planet.low_orbit);
        let (dv2, orbit) = orbit.with_circularization();

        //let dv0 = orbit.velocity_at_altitude(orbit.periapsis);

        let surface_velocity = (planet.radius * std::f64::consts::PI * 2.0) / planet.sidereal_day;
        let alt1 = planet.radius;
        let alt2 = planet.radius + planet.low_orbit;
        let gravity_loss = (2.0 * G * planet.mass * (alt2 - alt1) / (alt2 * alt1)).sqrt();

        //println!("{} {} {} {}", dv0, dv1, dv2, surface_velocity);

        (
            dv0 + dv1 + dv2 - surface_velocity + gravity_loss * 0.75,
            orbit,
        )
    }

    let mut orbit;
    let dv;

    (dv, orbit) = launch(&sun);
    println!("Sun surface -> 10,000km orbit: {dv:.1}");

    //let dv;
    //(dv, orbit) = orbit.with_apoapsis(earth.orbit.unwrap().1 - sun.radius);
    //println!("sun 10,000km -> earth intercept: {dv:.1}");

    let orbit_sun = orbit;

    let dv;
    (dv, orbit) = orbit.with_hohmann(earth.orbit.unwrap().1 - sun.radius);
    println!("sun 10,000km -> earth: {dv:.1}");

    println!("{}", calculate_escape_dv(&earth, 5_000.0));

    let dv;
    (dv, orbit) = orbit.with_hohmann(mars.orbit.unwrap().1 - sun.radius);
    println!("sun earth -> mars: {dv:.1}");

    let dv;
    (dv, orbit) = orbit.with_hohmann(jupiter.orbit.unwrap().1 - sun.radius);
    println!("sun mars -> jupiter: {dv:.1}");

    orbit = orbit_sun;

    let dv;
    (dv, orbit) = orbit.with_hohmann(earth.orbit.unwrap().1 - sun.radius);
    println!("sun 10,000km -> earth: {dv:.1}");

    let dv;
    (dv, orbit) = orbit.with_hohmann(jupiter.orbit.unwrap().1 - sun.radius);
    println!("sun earth -> jupiter: {dv:.1}");

    (_, orbit) = launch(&earth);
    let dv;
    (dv, orbit) = orbit.with_escape();
    println!("leo -> earth escape: {dv:.1}");

    let mut fake_delta_vs = vec![0.0; semimajors.len()];
    /*let mut fake_delta_vs =
    vec![
        0.0,
        302.82354760465154,
        132.249804747023,
        654.2987664479141,
        676.6456252145225,
        209.93262962578743,
        104.86155430590283,
        0.0,
    ];
        */
    let mut hops = vec![];

    for (i, a) in semimajors.iter().enumerate() {
        let mut orbit;
        let (dv_launch, dv_escape);
        (dv_launch, orbit) = launch(a.2.unwrap());

        (dv_escape, orbit) = orbit.with_escape();
        println!(
            "{: >10} {: >6.1}m/s launch Δv, {: >6.1}m/s escape Δv",
            a.0, dv_launch, dv_escape
        );
    }

    for (i, a) in semimajors.iter().enumerate() {
        let mut matrix_from_here = vec![];
        for b in &semimajors[i + 1..] {
            if a.1 == b.1 {
                continue;
            }

            let mut orbit;
            (_, orbit) = launch(&sun);
            (_, orbit) = orbit.with_hohmann(a.1 - sun.radius);
            let (mut dv1, mut dv2, orbit) = orbit.with_hohmann_2(b.1 - sun.radius);

            if let Some(planet) = a.2 {
                dv1 = calculate_escape_dv(planet, dv1);
            }

            if let Some(planet) = b.2 {
                dv2 = calculate_escape_dv(planet, dv2);
            }

            let dv = dv1 + dv2;

            println!(
                "{} -> {}: {:.1}m/s Δv ({:.1} + {:.1})",
                a.0, b.0, dv, dv1, dv2
            );

            matrix_from_here.push(dv);
        }

        hops.push((a.0, matrix_from_here));
    }

    fn dv_between_ground_truth(hops: &Hops, a: usize, b: usize) -> f64 {
        let (a, b) = (a.min(b), a.max(b));
        if b == a {
            return 0.0;
        }

        hops[a].1[b - a - 1]
    }

    fn show(
        fake_delta_vs: &[f64],
        hops: &Hops,
        semimajors: &Semimajors,
        print: bool,
    ) -> Option<f64> {
        if fake_delta_vs.iter().skip(1).any(|&x| x < 0.0) {
            //return None;
        }

        let mut score = 0.0;
        for (i, a) in semimajors.iter().enumerate() {
            for (j, b) in semimajors.iter().enumerate().skip(i + 1) {
                let mut hopping_dv = 0.0;
                for k in i..j {
                    hopping_dv += dv_between_ground_truth(hops, k, k + 1);
                    if k > i {
                        hopping_dv += fake_delta_vs[k];
                    }
                }

                hopping_dv *= 1.0 + fake_delta_vs[0] * 1.0;

                let ground_truth = dv_between_ground_truth(hops, i, j);
                let mut penalty = ((hopping_dv - ground_truth) / ground_truth).abs().powf(2.0);

                if hopping_dv < ground_truth {
                    penalty *= 1.0;
                }

                if i == 2 || j == 2 {
                    penalty = penalty.powf(3.0);
                    penalty *= 2000000000.0;
                }

                score += penalty;

                /*
                if hopping_dv < ground_truth {
                    return None;
                }
                */

                if print {
                    if i == 2 || j == 2 {
                        print!("\x1b[01;31m");
                    }
                    println!(
                    "{: >20}    {: >8.1} m/s Δv      {: >8.1} m/s Δv         {:>4.0}%   {:#>bar$}",
                    format!("{} → {}", a.0, b.0),
                    hopping_dv,
                    ground_truth,
                    (hopping_dv - ground_truth) / ground_truth * 100.0,
                    "",
                    bar = ((hopping_dv - ground_truth) / ground_truth * 100.0)
                        .floor()
                        .abs() as usize
                );
                    print!("\x1b[00m");
                }
            }
        }

        if print {
            for (i, (name, _, planet)) in semimajors.iter().enumerate() {
                print!("  {}", name);
                let scale = 1.0 + fake_delta_vs[0];
                if i < semimajors.len() - 1 {
                    if i > 0 && fake_delta_vs[i] != 0.0 {
                        print!("{:+.1}", fake_delta_vs[i] * scale);
                    }

                    let dv = dv_between_ground_truth(hops, i, i + 1);
                    print!(" {:.1}", dv * scale);
                }
            }
            println!();
        }

        Some(score)
    }

    let mut change = 1000.0;
    let mut index = 0;
    let mut any_worked = false;
    let mut iters: i32 = 5;
    loop {
        if iters >= 0 {
            iters -= 1;
        }
        if iters == 0 {
            println!("{:#?}", fake_delta_vs);
            show(&fake_delta_vs, &hops, semimajors, true);
        }

        let prev_score = show(&fake_delta_vs, &hops, semimajors, false).unwrap();
        fake_delta_vs[index] += change;
        let works = true;

        if let Some(score) = show(&fake_delta_vs, &hops, semimajors, false)
            && score < prev_score
        {
            if iters <= 0 {
                iters = 30000;
            }
            any_worked = true;
            continue;
        };

        fake_delta_vs[index] -= change;

        index += 1;
        if index >= semimajors.len() {
            index = 0;
            if any_worked {
                change *= 2.0;
            } else if change > 0.0 {
                change = -change;
            } else {
                change = -change * 0.8;
            }
            any_worked = false;
        }
    }

    println!(
        "LEO -> GTO ΔV: {:.1}",
        orbit_change_delta_v(
            &earth,
            250_000.0,
            250_000.0 + earth.radius,
            (250_000.0 + 35_768_000.0) / 2.0 + earth.radius,
        )
    );

    println!(
        "GTO -> GEO ΔV: {:.1}",
        orbit_change_delta_v(
            &earth,
            35_768_000.0,
            (250_000.0 + 35_768_000.0) / 2.0 + earth.radius,
            35_768_000.0 + earth.radius,
        )
    );
}
