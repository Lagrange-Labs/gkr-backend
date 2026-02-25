use crate::{Error, PCSFriParam, PolynomialCommitmentScheme, SecurityLevel};
use ff_ext::ExtensionField;
use multilinear_extensions::mle::ArcMultilinearExtension;
use witness::RowMajorMatrix;
use transcript::Transcript;


#[derive(Clone)]
pub struct MockPcs<F> {
    _marker: std::marker::PhantomData<F>,
}

impl PCSFriParam for () {
    fn get_pow_bits_by_level(&self, _pow_strategy: crate::PowStrategy) -> usize {
        1
    }
}
impl<E: ExtensionField> PolynomialCommitmentScheme<E> for MockPcs<E> {
    type Param = ();
    type ProverParam = ();
    type VerifierParam = ();
    type Commitment = ();
    type CommitmentChunk = ();
    type CommitmentWithWitness = ();
    type Proof = ();

    fn setup(_poly_size: usize, _security_level: SecurityLevel) -> Result<Self::Param, Error> {
        Ok(())
    }

    fn trim(
        _param: Self::Param,
        _poly_size: usize,
     ) -> Result<(Self::ProverParam, Self::VerifierParam), Error> {
        Ok(((), ()))
     }


    fn commit(
        _pp: &Self::ProverParam,
        _rmm: RowMajorMatrix<E::BaseField>,
     ) -> Result<Self::CommitmentWithWitness, Error> {
        Ok(())
     }

    fn batch_commit(
        _pp: &Self::ProverParam,
        _rmms: Vec<RowMajorMatrix<E::BaseField>>,
     ) -> Result<Self::CommitmentWithWitness, Error> {
        Ok(())
     }

    fn open(
        _pp: &Self::ProverParam,
        _poly: &ArcMultilinearExtension<E>,
        _comm: &Self::CommitmentWithWitness,
        _point: &[E],
        _eval: &E,
        _transcript: &mut impl Transcript<E>,
     ) -> Result<Self::Proof, Error> {
        Ok(())
     }

    fn verify(
        _vp: &Self::VerifierParam,
        _comm: &Self::Commitment,
        _point: &[E],
        _eval: &E,
        _proof: &Self::Proof,
        _transcript: &mut impl Transcript<E>,
    ) -> Result<(), Error> {
         Ok(())
     }
    
    fn write_commitment(
        _comm: &Self::Commitment,
        _transcript: &mut impl Transcript<E>,
    ) -> Result<(), Error> {
        Ok(())
    }
    
    fn get_pure_commitment(_comm: &Self::CommitmentWithWitness) -> Self::Commitment {
        ()
    }
    
    fn batch_open(
        _pp: &Self::ProverParam,
        _rounds: Vec<(
            &Self::CommitmentWithWitness,
            // for each matrix open at one point
            Vec<(crate::Point<E>, Vec<E>)>,
        )>,
        _transcript: &mut impl Transcript<E>,
    ) -> Result<Self::Proof, Error> {
        Ok(())
    }

    fn simple_batch_open(
        _pp: &Self::ProverParam,
        _polys: &[multilinear_extensions::mle::ArcMultilinearExtension<E>],
        _comm: &Self::CommitmentWithWitness,
        _point: &[E],
        _evals: &[E],
        _transcript: &mut impl Transcript<E>,
    ) -> Result<Self::Proof, Error> {
        Ok(())
    }
    
    fn batch_verify(
        _vp: &Self::VerifierParam,
        _rounds: Vec<(
            Self::Commitment,
            // for each matrix:
            Vec<(
                // its num_vars,
                usize,
                (
                    // the point,
                    crate::Point<E>,
                    // values at the point
                    Vec<E>,
                ),
            )>,
        )>,
        _proof: &Self::Proof,
        _transcript: &mut impl Transcript<E>,
    ) -> Result<(), Error> {
        Ok(())
    }
    
    fn simple_batch_verify(
        _vp: &Self::VerifierParam,
        _comm: &Self::Commitment,
        _point: &[E],
        _evals: &[E],
        _proof: &Self::Proof,
        _transcript: &mut impl Transcript<E>,
    ) -> Result<(), Error> {
        Ok(())
    }
    
    fn get_arc_mle_witness_from_commitment(
        _commitment: &Self::CommitmentWithWitness,
    ) -> Vec<multilinear_extensions::mle::ArcMultilinearExtension<'static, E>> {
        vec![]
    }
}

#[cfg(test)]
mod tests {
    use ff_ext::GoldilocksExt2 as F;
    use p3::field::FieldAlgebra;
    use rand::rngs::OsRng;
    use transcript::BasicTranscript;

    use super::*;

    #[test]
    fn test_mock_pcs() {
        let _ = MockPcs::<F>::setup(10, SecurityLevel::Conjecture100bits).unwrap();
        let (pp, vp) = MockPcs::<F>::trim((), 10).unwrap();
        let num_vars = 8;
        let rmm = RowMajorMatrix::<<F as ff_ext::ExtensionField>::BaseField>::rand(&mut OsRng, 1 << num_vars, 1);
        let poly: ArcMultilinearExtension<F> = rmm.to_mles().remove(0).into();
        let comm_with_witness = MockPcs::<F>::commit(&pp, rmm.clone()).unwrap();
        let comm = MockPcs::<F>::get_pure_commitment(&comm_with_witness);
        let mut transcript = BasicTranscript::new(b"hello");
        MockPcs::<F>::write_commitment(&comm, &mut transcript).unwrap();
        let rmms = vec![rmm.clone()];
        let _batch_comm = MockPcs::<F>::batch_commit(&pp, rmms).unwrap();
        let eval = F::ZERO;
        let point = vec![F::ZERO; 4];
        let proof = MockPcs::<F>::open(&pp, &poly, &comm_with_witness, &point, &eval, &mut transcript).unwrap();
        MockPcs::<F>::verify(&vp, &comm, &point, &eval, &proof, &mut transcript).unwrap();
        let rounds = vec![(
            &comm_with_witness,
            vec![(point.clone(), vec![eval])],
        )];
        let _batch_proof = MockPcs::<F>::batch_open(&pp, rounds, &mut transcript).unwrap();
        let polys = vec![poly.clone()];
        let _simple_proof = MockPcs::<F>::simple_batch_open(&pp, &polys, &comm_with_witness, &point, &[eval], &mut transcript).unwrap();
        let batch_rounds = vec![(comm, vec![(4, (point.clone(), vec![eval]))])];
        let _ = MockPcs::<F>::batch_verify(&vp, batch_rounds, &_batch_proof, &mut transcript);
        let _ = MockPcs::<F>::simple_batch_verify(&vp, &comm, &point, &[eval], &_simple_proof, &mut transcript);
        let _ = MockPcs::<F>::get_arc_mle_witness_from_commitment(&comm_with_witness);
    }
}