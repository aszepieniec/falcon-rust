use crate::falcon;

pub type SecretKey = falcon::SecretKey<1024>;
pub type PublicKey = falcon::PublicKey<1024>;
pub type Signature = falcon::Signature<1024>;

pub fn keygen(seed: [u8; 32]) -> (SecretKey, PublicKey) {
    falcon::keygen(seed)
}

pub fn sign(msg: &[u8], sk: &SecretKey) -> Signature {
    falcon::sign(msg, sk)
}

pub fn verify(msg: &[u8], sig: &Signature, pk: &PublicKey) -> bool {
    falcon::verify(msg, sig, pk)
}
