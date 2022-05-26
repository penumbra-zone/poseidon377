#[allow(non_snake_case)]
#[cfg(test)]
mod tests {
    use poseidon_paramgen::{Alpha, InputParameters, RoundNumbers};
    use std::convert::TryFrom;

    use ark_ff::{fields::FpParameters, BigInteger768};

    use ark_ed_on_bls12_377::FqParameters as Fq377Parameters;
    use num_bigint::BigUint;

    /// Represents a row in Table 7-9 in Appendix G of the paper.
    #[allow(dead_code)]
    struct TableRow {
        // Prime
        p: BigInteger768,
        // Security margin
        M: usize,
        // Text size: N = n * t
        N: usize,
        // S-box size: n
        n: usize,
        // Number of S-boxes (t)
        t: usize,
        // Number of full rounds
        r_F: usize,
        // Number of partial rounds
        r_P: usize,
        // Cost
        cost: usize,
    }

    #[test]
    fn table_7_concrete_instances_alpha_3() {
        let alpha = Alpha::Exponent(3);
        let N = 1536;

        // Appendix G of the paper provides in Table 7 concrete instances of x^{3} Poseidon
        let table_7 = [
            [128, N, 768, 2, 8, 83, 99],
            [128, N, 384, 4, 8, 84, 116],
            [128, N, 256, 6, 8, 84, 132],
            [128, N, 192, 8, 8, 84, 148],
            [128, N, 96, 16, 8, 64, 192],
            [256, N, 768, 2, 8, 170, 186],
            [256, N, 384, 4, 8, 171, 203],
            [256, N, 256, 6, 8, 171, 219],
            [256, N, 192, 8, 8, 128, 192],
            [256, N, 96, 16, 8, 64, 192],
        ];

        // These were randomly generated using `Crypto.Util.number.getPrime` in `calc_round_numbers.py`.
        // Note that changing these primes can shift the round numbers by one due to the fact that
        // in the Rust code we minimize `(M, log2(p))` instead of `(M, n)` for the interpolation and Grobner
        // constraints - as the equations in the text do. The code used to generate the LaTeX tables in the paper
        // minimized `(M, n)` where `n=ceil(log2(p))`.
        let primes = [
            "1285972319676842810969485878596785460378254587900320538194549286320222523621753054582270692761423067723622289509246331465441790711687096150336033426584236789046117244514940542370278948644966638168494119742350240132654417197013898277",
            "24029335548108648199755136672461166148126852491413827973546126772092090492121549469372885934460985132012930491106803",
            "92378444953952135362233764920562443588943051860418064632674652729586675198499",
            "5493676949319043689334439734705363934771468973473846167831",
            "47080328286681086157809746403",
            "1429028743336383885086224938533526527537492156009800001314231760304930487145108908025827879389380438318762771087581080106872120174685736157416469364705069428639169055604839082483182324797122682938663643974667493765095142241383410523",
            "39087297669469656432016004003288613036194633791610788101771094240290797024531975382555581671704748275265898438753551",
            "66599197946191249631812188829622192768573341198092727662820540460734999752353",
            "5939690934428747319089149475890099179736888602105582009677",
            "51558334253054077805045591617",
        ];

        for (row, prime) in table_7.iter().zip(primes.iter()) {
            let table_row = TableRow {
                p: BigInteger768::try_from(BigUint::parse_bytes(prime.as_bytes(), 10).unwrap())
                    .unwrap(),
                M: row[0],
                N: row[1],
                n: row[2],
                t: row[3],
                r_F: row[4],
                r_P: row[5],
                cost: row[6],
            };
            let input = InputParameters::new(table_row.M, table_row.t, table_row.p);
            let rounds = RoundNumbers::new(&input, &alpha);
            assert_eq!(rounds.full(), table_row.r_F);
            assert_eq!(rounds.partial(), table_row.r_P);
        }
    }

    #[test]
    fn table_8_concrete_instances_alpha_5() {
        let alpha = Alpha::Exponent(5);
        let N = 1536;

        // Appendix G of the paper provides in Table 8 concrete instances of x^{5} Poseidon
        let table_8 = [
            [128, N, 768, 2, 8, 56, 72],
            [128, N, 384, 4, 8, 56, 88],
            [128, N, 256, 6, 8, 57, 105],
            [128, N, 192, 8, 8, 57, 121],
            [128, N, 96, 16, 8, 42, 170],
            [256, N, 768, 2, 8, 116, 132],
            [256, N, 384, 4, 8, 116, 148],
            [256, N, 256, 6, 8, 117, 165],
            [256, N, 192, 8, 8, 86, 150],
            [256, N, 96, 16, 8, 42, 170],
        ];
        // These are randomly generated using `Crypto.Util.number.getPrime` in `calc_round_numbers.py`.
        let primes = [
            "1020432616714722447377130645995550938177200562461585225803048090516069775997125730278350623689059365786930362246398753740661742830734528163817673204838193731449505575391236215635076415262771576952755202348116656152449348660449600867",
            "39368751381006269902108472354351555399786464903134915790997916210058217900087282198897117382974492349917393640786091",
            "87461477668777205134073563425227256275751479514248606995724082059888688934739",
            "4848002283653915458342314818941918729698710728189008381529",
            "75652816395934777081987765559",
            "856194484071451972395636990297141234782681385983937272053939087749607049617620008843035758810491502978646254234829526865512106203220368169398938929573266422733785714524574286913159475045404392812568690255695103992898509688902585129",
            "31264640049088477770225752361009202019421489383518901149815764179913676105652492722523545538621601951766810698065469",
            "105572369273569957500357375582522416170139120006790850453866654215794611603289",
            "4158596063100046037384281156248765277093903050809579730379",
            "69186314987059605426149080777",
        ];

        for (row, prime) in table_8.iter().zip(primes.iter()) {
            let table_row = TableRow {
                p: BigInteger768::try_from(BigUint::parse_bytes(prime.as_bytes(), 10).unwrap())
                    .unwrap(),
                M: row[0],
                N: row[1],
                n: row[2],
                t: row[3],
                r_F: row[4],
                r_P: row[5],
                cost: row[6],
            };

            let input = InputParameters::new(table_row.M, table_row.t, table_row.p);
            let rounds = RoundNumbers::new(&input, &alpha);
            assert_eq!(rounds.full(), table_row.r_F);
            assert_eq!(rounds.partial(), table_row.r_P);
        }
    }

    #[test]
    fn table_9_concrete_instances_inverse_alpha() {
        let alpha = Alpha::Inverse;
        let N = 1536;

        // Appendix G of the paper provides in Table 9 concrete instances of x^{-1} Poseidon
        let table_9 = [
            [128, N, 768, 2, 8, 65, 81],
            [128, N, 384, 4, 8, 60, 92],
            [128, N, 256, 6, 8, 57, 105],
            [128, N, 192, 8, 8, 54, 118],
            [128, N, 96, 16, 8, 32, 160],
            [256, N, 768, 2, 8, 134, 150],
            [256, N, 384, 4, 8, 128, 160],
            [256, N, 256, 6, 8, 126, 174],
            [256, N, 192, 8, 8, 89, 153],
            [256, N, 96, 16, 8, 32, 160],
        ];
        // These are randomly generated using `Crypto.Util.number.getPrime` in `calc_round_numbers.py`.
        let primes = [
            "839024110632800606005488589415525259123742723495346937519809084766265078197349614971418834379058183972802414021198755701927311304539635955259602182310638445617111729043988505713791611514175087002304352201298897920438489996991375647",
            "36696523052796762519779528815436247683733002304013097671753885192987785136095786761462507148976873850586188951814571",
            "75170472159102037257875433833623287002137189164830086604962472127056704648801",
            "4894876304672917647761176180800568630532679564002326051649",
            "55291203714785298131871609871",
            "927284578876081672747839788289663410827541025085878755210392958787486165905270864807043744651809574495344863696391615387707402468149752052529827153626345801516239736053994574559414956937455670272070975463891435366211380843313541253",
            "31760045684781156881012718708292985274641612979396951420346125470291021841243802748526887313612061222985166447340549",
            "107799350998544074451708570201428080568984443711456541988628636420463500656153",
            "5784493250534137918494041902636593642197528046951094907327",
            "58525606832309595838417113617"
        ];

        for (row, prime) in table_9.iter().zip(primes.iter()) {
            let table_row = TableRow {
                p: BigInteger768::try_from(BigUint::parse_bytes(prime.as_bytes(), 10).unwrap())
                    .unwrap(),
                M: row[0],
                N: row[1],
                n: row[2],
                t: row[3],
                r_F: row[4],
                r_P: row[5],
                cost: row[6],
            };
            let input = InputParameters::new(table_row.M, table_row.t, table_row.p);
            let rounds = RoundNumbers::new(&input, &alpha);
            assert_eq!(rounds.full(), table_row.r_F);
            assert_eq!(rounds.partial(), table_row.r_P);
        }
    }

    #[test]
    fn poseidon_bls12_377_instance() {
        let alpha = Alpha::Exponent(17);

        // $t=3$ corresponds to a 2:1 hash
        let input = InputParameters::new(128, 3, Fq377Parameters::MODULUS);
        let rounds = RoundNumbers::new(&input, &alpha);
        assert_eq!(rounds.full(), 8);
        assert_eq!(rounds.partial(), 31);

        // $t=4$ corresponds to a 3:1 hash
        let input = InputParameters::new(128, 4, Fq377Parameters::MODULUS);
        let rounds = RoundNumbers::new(&input, &alpha);
        assert_eq!(rounds.full(), 8);
        assert_eq!(rounds.partial(), 31);

        // $t=5$ corresponds to a 4:1 hash
        let input = InputParameters::new(128, 5, Fq377Parameters::MODULUS);
        let rounds = RoundNumbers::new(&input, &alpha);
        assert_eq!(rounds.full(), 8);
        assert_eq!(rounds.partial(), 31);

        // $t=6$ corresponds to a 5:1 hash
        let input = InputParameters::new(128, 6, Fq377Parameters::MODULUS);
        let rounds = RoundNumbers::new(&input, &alpha);
        assert_eq!(rounds.full(), 8);
        assert_eq!(rounds.partial(), 31);
    }
}
