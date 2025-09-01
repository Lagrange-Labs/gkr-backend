pub mod impl_goldilocks {
    use std::sync::LazyLock;

    use crate::{
        ExtensionField, FieldFrom, FieldInto, FromUniformBytes, SmallField,
        array_try_from_uniform_bytes, impl_from_uniform_bytes_for_binomial_extension,
        poseidon::PoseidonField,
    };
    #[cfg(debug_assertions)]
    use p3::goldilocks::Poseidon2GoldilocksHL;
    use p3::{
        challenger::DuplexChallenger,
        field::{
            Field, FieldAlgebra, FieldExtensionAlgebra, PackedValue, PrimeField64, TwoAdicField,
            extension::{BinomialExtensionField, BinomiallyExtendable},
        },
        goldilocks::{
            Goldilocks, HL_GOLDILOCKS_8_EXTERNAL_ROUND_CONSTANTS,
            HL_GOLDILOCKS_8_INTERNAL_ROUND_CONSTANTS, MATRIX_DIAG_8_GOLDILOCKS,
        },
        merkle_tree::MerkleTreeMmcs,
        poseidon2::{
            ExternalLayer, InternalLayer, MDSMat4, add_rc_and_sbox_generic,
            external_initial_permute_state, external_terminal_permute_state,
            internal_permute_state, matmul_internal,
        },
        symmetric::{PaddingFreeSponge, Permutation, TruncatedPermutation},
    };

    #[cfg(debug_assertions)]
    use crate::poseidon::impl_instruments::*;

    use p3::symmetric::CryptographicPermutation;

    pub type GoldilocksExt2 = BinomialExtensionField<Goldilocks, 2>;

    impl FieldFrom<u64> for Goldilocks {
        fn from_v(v: u64) -> Self {
            Self::from_canonical_u64(v)
        }
    }

    impl FieldFrom<u64> for GoldilocksExt2 {
        fn from_v(v: u64) -> Self {
            Self::from_canonical_u64(v)
        }
    }

    impl FieldInto<Goldilocks> for Goldilocks {
        fn into_f(self) -> Goldilocks {
            self
        }
    }

    pub const POSEIDON2_GOLDILICK_WIDTH: usize = 8;
    pub const POSEIDON2_GOLDILICK_RATE: usize = 4;

    #[cfg(debug_assertions)]
    impl CryptographicPermutation<[Goldilocks; POSEIDON2_GOLDILICK_WIDTH]>
        for Instrumented<Poseidon2GoldilocksHL<POSEIDON2_GOLDILICK_WIDTH>>
    {
    }
    /// Implements the poseidon permutation without the need for allocations.
    #[derive(Copy, Clone)]
    pub struct NoAllocPoseidon {}

    #[derive(Copy, Clone)]
    struct NoAllocExternalLayer();

    #[derive(Copy, Clone)]
    struct NoAllocInternalLayer();

    const WIDTH: usize = 8;
    const GOLDILOCKS_S_BOX_DEGREE: u64 = 7;
    static INITIAL_EXTERNAL_CONSTANTS: LazyLock<[[Goldilocks; WIDTH]; 4]> = LazyLock::new(|| {
        HL_GOLDILOCKS_8_EXTERNAL_ROUND_CONSTANTS[0]
            .map(|inner| inner.map(Goldilocks::from_canonical_u64))
    });
    static TERMINAL_EXTERNAL_CONSTANTS: LazyLock<[[Goldilocks; WIDTH]; 4]> = LazyLock::new(|| {
        HL_GOLDILOCKS_8_EXTERNAL_ROUND_CONSTANTS[1]
            .map(|inner| inner.map(Goldilocks::from_canonical_u64))
    });
    static INTERNAL_CONSTANTS: LazyLock<[Goldilocks; 22]> = LazyLock::new(|| {
        HL_GOLDILOCKS_8_INTERNAL_ROUND_CONSTANTS.map(Goldilocks::from_canonical_u64)
    });

    impl<FA: FieldAlgebra<F = Goldilocks>> ExternalLayer<FA, WIDTH, GOLDILOCKS_S_BOX_DEGREE>
        for NoAllocExternalLayer
    {
        fn permute_state_initial(&self, state: &mut [FA; WIDTH]) {
            external_initial_permute_state(
                state,
                &*INITIAL_EXTERNAL_CONSTANTS,
                add_rc_and_sbox_generic::<_, GOLDILOCKS_S_BOX_DEGREE>,
                &MDSMat4,
            );
        }

        fn permute_state_terminal(&self, state: &mut [FA; WIDTH]) {
            external_terminal_permute_state(
                state,
                &*TERMINAL_EXTERNAL_CONSTANTS,
                add_rc_and_sbox_generic::<_, GOLDILOCKS_S_BOX_DEGREE>,
                &MDSMat4,
            );
        }
    }

    impl<FA: FieldAlgebra<F = Goldilocks>> InternalLayer<FA, WIDTH, GOLDILOCKS_S_BOX_DEGREE>
        for NoAllocInternalLayer
    {
        /// Perform the internal layers of the Poseidon2 permutation on the given state.
        fn permute_state(&self, state: &mut [FA; 8]) {
            internal_permute_state::<FA, 8, GOLDILOCKS_S_BOX_DEGREE>(
                state,
                |x| matmul_internal(x, MATRIX_DIAG_8_GOLDILOCKS),
                &*INTERNAL_CONSTANTS,
            )
        }
    }

    static EXTERNAL_LAYER: NoAllocExternalLayer = NoAllocExternalLayer();
    static INTERNAL_LAYER: NoAllocInternalLayer = NoAllocInternalLayer();

    impl Permutation<[Goldilocks; WIDTH]> for NoAllocPoseidon {
        fn permute_mut(&self, state: &mut [Goldilocks; WIDTH]) {
            EXTERNAL_LAYER.permute_state_initial(state);
            INTERNAL_LAYER.permute_state(state);
            EXTERNAL_LAYER.permute_state_terminal(state);
        }
    }

    impl CryptographicPermutation<[Goldilocks; WIDTH]> for NoAllocPoseidon {}

    #[cfg(debug_assertions)]
    impl CryptographicPermutation<[Goldilocks; WIDTH]> for Instrumented<NoAllocPoseidon> {}

    impl PoseidonField for Goldilocks {
        #[cfg(debug_assertions)]
        type P = Instrumented<NoAllocPoseidon>;
        #[cfg(not(debug_assertions))]
        type P = NoAllocPoseidon;
        type T =
            DuplexChallenger<Self, Self::P, POSEIDON2_GOLDILICK_WIDTH, POSEIDON2_GOLDILICK_RATE>;
        type S = PaddingFreeSponge<Self::P, POSEIDON2_GOLDILICK_WIDTH, POSEIDON2_GOLDILICK_RATE, 4>;
        type C = TruncatedPermutation<Self::P, 2, 4, POSEIDON2_GOLDILICK_WIDTH>;
        type MMCS = MerkleTreeMmcs<Self, Self, Self::S, Self::C, 4>;
        fn get_default_challenger() -> Self::T {
            DuplexChallenger::<Self, Self::P, POSEIDON2_GOLDILICK_WIDTH, POSEIDON2_GOLDILICK_RATE>::new(
                Self::get_default_perm(),
            )
        }

        #[cfg(debug_assertions)]
        fn get_default_perm() -> Self::P {
            Instrumented::new(NoAllocPoseidon {})
        }

        #[cfg(not(debug_assertions))]
        fn get_default_perm() -> Self::P {
            NoAllocPoseidon {}
        }

        fn get_default_sponge() -> Self::S {
            PaddingFreeSponge::new(Self::get_default_perm())
        }

        fn get_default_compression() -> Self::C {
            TruncatedPermutation::new(Self::get_default_perm())
        }

        fn get_default_mmcs() -> Self::MMCS {
            MerkleTreeMmcs::new(Self::get_default_sponge(), Self::get_default_compression())
        }
    }

    impl_from_uniform_bytes_for_binomial_extension!(p3::goldilocks::Goldilocks, 2);

    impl FromUniformBytes for Goldilocks {
        type Bytes = [u8; 8];

        fn try_from_uniform_bytes(bytes: [u8; 8]) -> Option<Self> {
            let value = u64::from_le_bytes(bytes);
            let is_canonical = value < Self::ORDER_U64;
            is_canonical.then(|| Self::from_canonical_u64(value))
        }
    }

    impl SmallField for Goldilocks {
        const MODULUS_U64: u64 = Self::ORDER_U64;

        /// Convert a byte string into a list of field elements
        fn bytes_to_field_elements(bytes: &[u8]) -> Vec<Self> {
            bytes
                .chunks(8)
                .map(|chunk| {
                    let mut array = [0u8; 8];
                    array[..chunk.len()].copy_from_slice(chunk);
                    unsafe { std::ptr::read_unaligned(array.as_ptr() as *const u64) }
                })
                .map(Self::from_canonical_u64)
                .collect::<Vec<_>>()
        }

        /// Convert a field elements to a u64.
        fn to_canonical_u64(&self) -> u64 {
            self.as_canonical_u64()
        }
    }

    impl ExtensionField for GoldilocksExt2 {
        const DEGREE: usize = 2;
        const MULTIPLICATIVE_GENERATOR: Self = <GoldilocksExt2 as Field>::GENERATOR;
        const TWO_ADICITY: usize = Goldilocks::TWO_ADICITY;
        // non-residue is the value w such that the extension field is
        // F[X]/(X^2 - w)
        const NONRESIDUE: Self::BaseField = <Goldilocks as BinomiallyExtendable<2>>::W;

        type BaseField = Goldilocks;

        fn to_canonical_u64_vec(&self) -> Vec<u64> {
            self.as_base_slice()
                .iter()
                .map(|v: &Self::BaseField| v.as_canonical_u64())
                .collect()
        }
    }
}
