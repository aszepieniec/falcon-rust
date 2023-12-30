use std::{
    fs::File,
    io::{Read, Write},
};

use falcon_rust::falcon::{PublicKey, SecretKey, Signature, FALCON_1024, FALCON_512};
use rand::{thread_rng, Rng};

pub fn save_helper_data() -> std::io::Result<()> {
    let mut rng = thread_rng();
    let (sk512, pk512) = FALCON_512.keygen(rng.gen());
    let msg = "Hello, world!".as_bytes();
    let sig512 = FALCON_512.sign(msg, &sk512);

    File::create("helper_data/sk512.dat")?.write_all(&sk512.to_bytes())?;
    File::create("helper_data/pk512.dat")?.write_all(&pk512.to_bytes())?;
    File::create("helper_data/msg512.dat")?.write_all(msg)?;
    File::create("helper_data/sig512.dat")?.write_all(&sig512.to_bytes())?;

    let mut rng = thread_rng();
    let (sk1024, pk1024) = FALCON_512.keygen(rng.gen());
    let msg = "Hello, world!".as_bytes();
    let sig1024 = FALCON_1024.sign(msg, &sk1024);

    File::create("helper_data/sk1024.dat")?.write_all(&sk1024.to_bytes())?;
    File::create("helper_data/pk1024.dat")?.write_all(&pk1024.to_bytes())?;
    File::create("helper_data/msg1024.dat")?.write_all(msg)?;
    File::create("helper_data/sig1024.dat")?.write_all(&sig1024.to_bytes())?;

    Ok(())
}

#[allow(clippy::type_complexity)]
pub fn load_helper_data() -> std::io::Result<(
    SecretKey,
    PublicKey,
    Vec<u8>,
    Signature,
    SecretKey,
    PublicKey,
    Vec<u8>,
    Signature,
)> {
    let mut buf = vec![];
    File::open("helper_data/sk512.dat")?.read_to_end(&mut buf)?;
    let sk512 = SecretKey::from_bytes(&buf).expect("Could not parse secret key.");
    File::open("helper_data/pk512.dat")?.read_to_end(&mut buf)?;
    let pk512 = PublicKey::from_bytes(&buf).expect("Could not parse public key.");
    File::open("helper_data/msg512.dat")?.read_to_end(&mut buf)?;
    let msg512 = buf.clone();
    File::open("helper_data/sig512.dat")?.read_to_end(&mut buf)?;
    let sig512 = Signature::from_bytes(&buf).expect("Could not parse signature.");

    File::open("helper_data/sk1024.dat")?.read_to_end(&mut buf)?;
    let sk1024 = SecretKey::from_bytes(&buf).expect("Could not parse secret key.");
    File::open("helper_data/pk1024.dat")?.read_to_end(&mut buf)?;
    let pk1024 = PublicKey::from_bytes(&buf).expect("Could not parse public key.");
    File::open("helper_data/msg1024.dat")?.read_to_end(&mut buf)?;
    let msg1024 = buf.clone();
    File::open("helper_data/sig1024.dat")?.read_to_end(&mut buf)?;
    let sig1024 = Signature::from_bytes(&buf).expect("Could not parse signature.");

    Ok((
        sk512, pk512, msg512, sig512, sk1024, pk1024, msg1024, sig1024,
    ))
}
