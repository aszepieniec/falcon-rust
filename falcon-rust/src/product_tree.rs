use lazy_static::lazy_static;
use std::collections::VecDeque;

use itertools::Itertools;
use num::BigUint;
use num::FromPrimitive;
use num::ToPrimitive;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProductBranch {
    pub product: BigUint,
    pub left: Box<ProductTree>,
    pub right: Box<ProductTree>,
}

impl ProductBranch {
    /// Generates a rust expression to construct an object identical to the
    /// given one. This function was used to generate the large `lazy_static`
    /// expression populating `MASTER_TREE`.
    #[allow(dead_code)]
    fn as_rust_expression(&self) -> String {
        format!(
            "ProductBranch {{ product: BigUint::from_slice(&[{}]), left: Box::new({}), right: Box::new({}) }} ",
            self.product.to_u32_digits().iter().join(", "),
            self.left.as_rust_expression(),
            self.right.as_rust_expression()
        )
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ProductTree {
    Leaf(u32),
    Branch(ProductBranch),
}

impl ProductTree {
    pub fn value(&self) -> BigUint {
        match self {
            ProductTree::Leaf(leaf) => BigUint::from_u32(*leaf).unwrap(),
            ProductTree::Branch(branch) => branch.product.clone(),
        }
    }

    /// Create a [`ProductTree`] from a list of leafs, which are primes fitting
    /// in u32s. This function was used to generate the large `lazy_static`
    /// expression populating `MASTER_TREE`.
    #[allow(dead_code)]
    pub fn from_leafs(integers: &[u32]) -> ProductTree {
        let mut deque: VecDeque<ProductTree> = integers
            .into_iter()
            .map(|&i| ProductTree::Leaf(i))
            .collect();
        while deque.len() >= 2 {
            let left = deque.pop_front().unwrap();
            let right = deque.pop_front().unwrap();
            let new_node = ProductTree::Branch(ProductBranch {
                product: left.value() * right.value(),
                left: Box::new(left),
                right: Box::new(right),
            });
            deque.push_back(new_node);
        }
        deque[0].clone()
    }

    pub fn reduce(&self, int: &BigUint) -> Vec<u32> {
        match self {
            ProductTree::Leaf(_leaf) => vec![int.to_u32().unwrap()],
            ProductTree::Branch(branch) => {
                let left_remainder = int % branch.left.value();
                let right_remainder = int % branch.right.value();
                [
                    branch.left.reduce(&left_remainder),
                    branch.right.reduce(&right_remainder),
                ]
                .concat()
            }
        }
    }

    pub fn as_branch(&self) -> &ProductBranch {
        match self {
            ProductTree::Leaf(_leaf) => {
                panic!("cannot cast product tree as branch when it is a leaf")
            }
            ProductTree::Branch(branch) => branch,
        }
    }

    pub fn as_rust_expression(&self) -> String {
        match self {
            ProductTree::Leaf(leaf) => format!("ProductTree::Leaf({})", *leaf),
            ProductTree::Branch(branch) => {
                format!("ProductTree::Branch({})", branch.as_rust_expression())
            }
        }
    }
}

lazy_static! {
    pub(crate) static ref MASTER_TREE: ProductTree = ProductTree::Branch(ProductBranch {
        product: BigUint::from_slice(&[
            3753484289, 1373989043, 2334866072, 3329235467, 3329256830, 1264736616, 1819217197,
            3209886709, 2256655030, 3612497141, 1391733266, 3945799240, 211919685, 3234670912,
            2627652267, 4083418222, 1323166918, 171708379, 2967545405, 2057090398, 1971047875,
            2128327237, 804477584, 905195573, 3507430947, 2417028454, 1979123497, 3387294551,
            2732020476, 1436551177, 897196786, 56157073, 1777587022, 3377542594, 133469762,
            1357171266, 68693498, 2945716239, 1995157799, 563452935, 3172771053, 4025754421,
            1375469350, 3365679907, 2880401243, 1672773333, 820804901, 1656038447, 719083955,
            1218826358, 3414817828, 2864370751, 582519587, 3582713367, 1059553150, 1135528286,
            1110474456, 1726455491, 2955744463, 888904401, 4036858668, 2412330920, 3598177228,
            2031327558, 766894109, 1155593834, 3927773671, 2971805235, 2427218523, 170908827,
            2825024315, 1639312291, 2368150649, 3145639915, 4244421236, 715693343, 2726834788,
            3863635653, 549461169, 2321494121, 2938327794, 105215796, 4177559220, 3060322617,
            1428393167, 2919550068, 2955145215, 2692395811, 3848183446, 755732362, 3747534829,
            4261624419, 710135847, 3265731089, 298787419, 865108383, 2662613763, 659895096,
            3937219409, 894199727, 2836695889, 1706338787, 2658608949, 570557051, 693775132,
            124674601, 2227774512, 1530818218, 3187138482, 1820594287, 3691956364, 186818138,
            367557946, 1308033737, 3987083143, 938984631, 3248634441, 1191233608, 3652674558,
            3967772334, 468108721, 1382490035, 254027183, 3085687493, 25421885, 3683205690,
            4194299491, 3651136695, 1479646240, 4059197836, 2917598002, 138606466, 2518711665,
            1100986715, 1224821022, 3463644102, 3855651540, 1100080399, 2708651155, 2845809953,
            2281316301, 1998474447, 3941305086, 1968408231, 2346799601, 3958519296, 1514785825,
            785409317, 1206049668, 3342017532, 2808091859, 732435221, 1790275930, 1944114102,
            490016712, 1623199045, 3984068061, 2116536175, 1440937961, 3806970387, 2340074014,
            2977348389, 3732964998, 1696603263, 3353934297, 1581552272, 1293711318, 4075462417,
            1610485794, 815658551, 3182230984, 232318934, 3356288495, 1830866737, 2922385473,
            2338372206, 2864112278, 3954262123, 2636002261, 2232429137, 3109725443, 3030383429,
            307626671, 4216078167, 3496543274, 1122560840, 3608619558, 1543281852, 1533163564,
            750200702, 2477366936, 2523981046, 172549400, 3431752597, 2403186452, 699722204,
            1125135316, 130119787, 3848800081, 874003430, 2097931707, 3989123972, 2860327920,
            2699739268, 3877075985, 3676041157, 3201169775, 1589194723, 2822556598, 2235948879,
            3036972173, 3293557853, 2616743509, 2643526932, 3265466485, 3849057931, 3107001766,
            3913548832, 1691835893, 3163192972, 1050408209, 2553408144, 1838035350, 3670902646,
            2539530794, 1540973897, 3677040464, 2807008548, 65058060, 3972297718, 2943219373,
            3011441265, 1937314911, 4138090352, 2153290941, 2009075152, 1509457009, 1109269361,
            74718770, 1717436524, 3885151250, 3470013757, 1049327802, 3610488143, 2976408870,
            1016988784, 1100572401, 2888581464, 520470801, 1942272256, 2061109335, 1381964062,
            1772652007, 3554459501, 2070668406, 3768749148, 1301659748, 3954874549, 3417395229,
            740657542, 2297577025, 2745122739, 803443450, 1219482499, 1341923288, 402317839,
            3522944421, 1676361980, 1401669063, 599280747, 2378407111, 3854073894, 3804410996,
            1655926216, 4020583923, 2765349054, 1389283608, 4075475019, 718210554, 2548662514,
            3359385583, 3182654289, 2494534957, 1588973760, 1341169812, 1482717319, 1178978913,
            523806030, 604297125, 4243392054, 2506354766, 3427683216, 1877251525, 3372145078,
            3234521570, 1953686076, 39454406, 769604162, 1261976333, 976173563, 3098928913,
            2346952447, 2042947550, 2747644871, 1602994885, 4153665401, 1351338767, 606418903,
            3571635102, 1343258016, 681706783, 2025919624, 1633505881, 2304059681, 1065669249,
            3185849431, 887796628, 1828127855, 347619059, 1294909748, 2139668914, 4186391518,
            3524401432, 4045154220, 3847302944, 2011736482, 336648795, 2447659100, 2700111442,
            2479374575, 3639727661, 1674314457, 4097854557, 4102319262, 3707163447, 1646133357,
            1658672087, 2962366492, 369094157, 3019397532, 3323623830, 500998056, 522672053,
            1092642211, 1878313983, 1663357236, 3197373258, 399842492, 3871349878, 4031932980,
            3988723017, 2050127857, 2535276770, 2558561082, 3054392927, 3063972230, 429526710,
            2740312473, 2982189680, 3691140868, 2541423701, 254113864, 2487443853, 3448727415,
            3278839931, 3357828444, 2493036160, 3763990232, 2388275359, 751925508, 3613791536,
            284275629, 3375090844, 4075698881, 2777897985, 2855026020, 2397879350, 2455292332,
            263275971, 1376211758, 876104029, 359706849, 564482311, 1622030107, 833019698,
            2261148756, 876482464, 3942074084, 3128834743, 3409870188, 3125782553, 3674250969,
            2252933514, 1422433444, 3059822248, 446691010, 798408011, 624440226, 1078221209,
            3646522598, 4264921157, 1045953362, 2154833470, 1538551799, 1106407303, 4024391734,
            1885028242, 735527206, 2292173339, 694101214, 2184482402, 445618626, 3486188702,
            2837196200, 2263752205, 3985319299, 3526183367, 588756941, 3721701059, 615979578,
            2952242738, 2509680601, 3790013072, 2551994794, 2104089399, 3082202168, 3947835627,
            1266421962, 1194041935, 1817560508, 3745808390, 1375740773, 1727471161, 4157879576,
            944101440, 1377350251, 4102026405, 3644630624, 51539423, 425552253, 2466831970,
            4137194539, 3413246538, 1815483460, 1652667693, 803089729, 2417410005, 2536359686,
            668324147, 1060094179, 3498635406, 1873361273, 745559152, 3373559529, 765377997,
            1597831835, 971277575, 41024555, 2814679293, 3116046217, 3207762830, 3965212282,
            242453576, 658503659, 2087307694, 1406162818, 1599073581, 3182305997, 163340118,
            2054372578, 1484427403, 4292413113, 3412934640, 1667579233, 3823961493, 3725957166,
            795983444, 3524432968, 3049088410, 720391504, 23879
        ]),
        left: Box::new(ProductTree::Branch(ProductBranch {
            product: BigUint::from_slice(&[
                9101313, 4121361246, 3438653410, 2384310242, 1745903058, 3372910688, 1553785176,
                714675419, 1369034381, 878936568, 3398144615, 1692840865, 3356171854, 2016584814,
                574267290, 3967267273, 4289121288, 784992127, 2653625434, 3799347477, 894191061,
                4185566945, 1135833169, 183263887, 1560591179, 1457615876, 4109472505, 579624628,
                3158211903, 3148903179, 2565456482, 819797829, 2850455428, 2501177652, 2413920824,
                936514793, 1236538614, 3343985662, 1669998930, 903961213, 4187774301, 2373428999,
                219342235, 3184207590, 987667790, 3325141256, 1623274815, 586888936, 2321502562,
                818606316, 2025095897, 1082370671, 70561051, 1140299787, 3307050099, 3800502492,
                4082618473, 3590805657, 1187622359, 3642765476, 255983270, 24210953, 3395233872,
                43064457, 2921090479, 3436876451, 2643210006, 2147598304, 4288103251, 280164141,
                1555752917, 2299096378, 4115022535, 915798349, 128056218, 2013991973, 1514217749,
                2010876600, 224605871, 1293296234, 1999599775, 3774155679, 523855548, 1014330413,
                143233471, 206815068, 1270452650, 3219459255, 575427732, 1618818529, 1829845188,
                1057332026, 1659429312, 3090028452, 3224985307, 3142577655, 1269020362, 2743784803,
                3915898875, 474269684, 752917847, 3515464271, 1972538186, 3334362038, 957784879,
                3251926655, 749581758, 1336586184, 3686877922, 210287002, 2937530328, 3405103412,
                367206166, 2591498165, 2337273236, 16524949, 246624449, 1219155030, 1067696553,
                3788764984, 3844829412, 1208990906, 3308258626, 1814834906, 3253135636, 343506373,
                1650940159, 1213005132, 4820972, 65059664, 3808751727, 2516248137, 901509300,
                3072989702, 581002958, 106040114, 1849051906, 396724879, 196710905, 1479196197,
                3766191264, 1734966930, 1864217717, 2042049171, 2217490783, 412206435, 409203423,
                1836089168, 1300536359, 3724544900, 1028788657, 2427664629, 23166243, 1164580391,
                3515621549, 2426337385, 4087776447, 77606077, 2885479221, 3548417089, 653313597,
                660393473, 4293538382, 3084647657, 1436265436, 629726739, 804637413, 1224933060,
                3475549331, 332384737, 1586086691, 3792240888, 1558549229, 1964900140, 2360601499,
                3889156314, 82102661, 2175649904, 2212082145, 977727487, 1441333936, 2230493831,
                413171496, 3837248694, 3743432613, 4044273980, 1859254178, 2940057989, 2469113917,
                817903459, 843953017, 2104697941, 3581187803, 61565559, 2998714895, 1321157085,
                1547443668, 3110781802, 577220791, 240985692, 1243974532, 3335698698, 3821962986,
                2028260412, 2366559629, 3139491581, 32277371, 2893279965, 594459473, 337572640,
                779713284, 1948548796, 738124548, 2117400553, 3209873920, 731571996, 1500446470,
                2323834186, 3291630387, 544790278, 1049930914, 1540020396, 4145369821, 1507847466,
                3545892279, 3099771735, 1178968405, 2428928255, 439388119, 819484993, 1199864449,
                28009453, 21083092, 1324370674, 3483906145, 786647414, 2542470783, 3872735229,
                3806487816, 1985110701, 12
            ]),
            left: Box::new(ProductTree::Branch(ProductBranch {
                product: BigUint::from_slice(&[
                    2065170433, 1110546307, 799182345, 3267657279, 2208294592, 3015640450,
                    1614536673, 85293852, 842879171, 3694726100, 2911769173, 1376637800, 516868310,
                    3986238837, 1133866458, 1520674534, 2616177212, 1966445954, 3458221298,
                    1594062773, 1455017941, 2366775067, 3351335306, 3448382563, 876641374,
                    4239803228, 930770731, 591010958, 1612854845, 553066135, 1660409591,
                    3537605525, 1463820484, 3046897394, 3786445558, 1012574097, 968209191,
                    2984691390, 426841893, 1594560872, 687572006, 390125276, 1259544576,
                    2371292241, 4255605680, 355724053, 1474241489, 3020195217, 200779018,
                    1565952487, 468808503, 4164013479, 576011957, 1162152959, 2685164763,
                    453385368, 2974745212, 1219360385, 2251383087, 1136504753, 2943760600,
                    2067498717, 368309998, 511024163, 2043160219, 842453927, 537441145, 123588610,
                    4144428305, 2332956086, 1673102475, 1666971212, 4140275818, 4015066146,
                    3806077419, 986418740, 101154748, 1398379439, 1685340741, 434927111, 159144854,
                    300687844, 3024314159, 302040165, 3558701800, 2213282968, 1338755895, 31236707,
                    470139146, 545952771, 1521723353, 3080523035, 3122283181, 3695540554,
                    617098353, 2342510111, 2475621578, 335363011, 3367239170, 3377771918,
                    944698623, 1284117446, 3078333087, 3071319603, 2850329586, 3707770886,
                    181797120, 3992254462, 2174332764, 67376533, 2248826241, 3354846729,
                    1189519079, 30558537, 2754185344, 4171246084, 4287644546, 1052878359,
                    2011822779, 3600835844, 1
                ]),
                left: Box::new(ProductTree::Branch(ProductBranch {
                    product: BigUint::from_slice(&[
                        2979741697, 2773963522, 1706700418, 3285706742, 3711196872, 2454309704,
                        1174935369, 3469514092, 3294070321, 1945162104, 3268995703, 1172685485,
                        1625326593, 3004441793, 1581147641, 3747167795, 3855971222, 3324367010,
                        1344082510, 2130790755, 3504819523, 602041594, 3271704986, 2234140351,
                        53812018, 4043598433, 3933930022, 2154560661, 1012995942, 2476971967,
                        1404410139, 1021496682, 712989272, 331455882, 3514677056, 2582102741,
                        1210912334, 3685062148, 3737667025, 3072308831, 2843292437, 2540452937,
                        2903029087, 743637444, 989127819, 2652687749, 2050501450, 513247981,
                        2991719566, 3875436692, 3973007060, 3568807309, 78613972, 1023187022,
                        1263082498, 3109496313, 2950169029, 2063822830, 1493648213, 694378728, 1
                    ]),
                    left: Box::new(ProductTree::Branch(ProductBranch {
                        product: BigUint::from_slice(&[
                            675569665, 279111285, 3569956568, 3129227589, 587425014, 2854600023,
                            1826021053, 1209368299, 4199712590, 1218822528, 3660417837, 2809247233,
                            2986422343, 3708926639, 659303050, 4007026270, 3538046992, 2074295087,
                            2286417007, 495443072, 1795006356, 3947049599, 2971591958, 250909550,
                            3897777684, 3592081072, 1449733492, 149699158, 2810848660, 154755008,
                            1
                        ]),
                        left: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                2426126337, 607456000, 491338359, 1100840650, 532376086, 504605603,
                                944938171, 2928371180, 4240212544, 139510643, 20072607, 3729602310,
                                3506681717, 1493572594, 41008619, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    3123732481, 2245682344, 2316416653, 1643649388, 4038227279,
                                    1824504148, 3163932821, 65729
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1879728129, 3843305669, 1460489281, 16787842
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[419651585, 268490753]),
                                        left: Box::new(ProductTree::Leaf(1073754113)),
                                        right: Box::new(ProductTree::Leaf(1073950721))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3070689281, 268550156]),
                                        left: Box::new(ProductTree::Leaf(1073958913)),
                                        right: Box::new(ProductTree::Leaf(1073983489))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        170262529, 1730398873, 2481311985, 16816161
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1024466945, 268699712]),
                                        left: Box::new(ProductTree::Leaf(1074196481)),
                                        right: Box::new(ProductTree::Leaf(1074343937))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[286646273, 268793976]),
                                        left: Box::new(ProductTree::Leaf(1074442241)),
                                        right: Box::new(ProductTree::Leaf(1074475009))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    2456510465, 648091533, 2502775218, 3296437043, 3017708231,
                                    2368056284, 3143434677, 65966
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        4264624129, 2830468228, 1461025216, 16827448
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2199371777, 268824717]),
                                        left: Box::new(ProductTree::Leaf(1074515969)),
                                        right: Box::new(ProductTree::Leaf(1074524161))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3407429633, 268849311]),
                                        left: Box::new(ProductTree::Leaf(1074548737)),
                                        right: Box::new(ProductTree::Leaf(1074589697))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        876240897, 3566872732, 3478944357, 16837071
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[387678209, 268886205]),
                                        left: Box::new(ProductTree::Leaf(1074597889)),
                                        right: Box::new(ProductTree::Leaf(1074688001))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[220127233, 268941550]),
                                        left: Box::new(ProductTree::Leaf(1074696193)),
                                        right: Box::new(ProductTree::Leaf(1074810881))
                                    }))
                                }))
                            }))
                        })),
                        right: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                1068015617, 4086531741, 398461313, 4246555490, 2569150060,
                                196145208, 3641540067, 4250579946, 1350313913, 1021702860,
                                1678628000, 1556288478, 155881197, 1123930958, 112670603, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    715841537, 155751711, 1507106894, 434398283, 2864331940,
                                    19483713, 2651033004, 66222
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        4132077569, 1697984926, 143450068, 16853634
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[388169729, 269009202]),
                                        left: Box::new(ProductTree::Leaf(1074835457)),
                                        right: Box::new(ProductTree::Leaf(1074941953))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1596424193, 269083014]),
                                        left: Box::new(ProductTree::Leaf(1075007489)),
                                        right: Box::new(ProductTree::Leaf(1075064833))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1147166721, 2060296500, 2056895885, 16876121
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3039649793, 269179391]),
                                        left: Box::new(ProductTree::Leaf(1075105793)),
                                        right: Box::new(ProductTree::Leaf(1075351553))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1865613313, 269271690]),
                                        left: Box::new(ProductTree::Leaf(1075376129)),
                                        right: Box::new(ProductTree::Leaf(1075449857))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    2365440001, 701150238, 1567995801, 3631112061, 4137414211,
                                    2045187643, 3837137936, 66557
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1014210561, 1611956979, 2450195020, 16895929
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1161273345, 269347589]),
                                        left: Box::new(ProductTree::Leaf(1075507201)),
                                        right: Box::new(ProductTree::Leaf(1075621889))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2000420865, 269419396]),
                                        left: Box::new(ProductTree::Leaf(1075695617)),
                                        right: Box::new(ProductTree::Leaf(1075720193))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2626297857, 2985147848, 705180699, 16919103
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2067906561, 269513780]),
                                        left: Box::new(ProductTree::Leaf(1075752961)),
                                        right: Box::new(ProductTree::Leaf(1076039681))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2034786305, 269622557]),
                                        left: Box::new(ProductTree::Leaf(1076064257)),
                                        right: Box::new(ProductTree::Leaf(1076162561))
                                    }))
                                }))
                            }))
                        }))
                    })),
                    right: Box::new(ProductTree::Branch(ProductBranch {
                        product: BigUint::from_slice(&[
                            1968627713, 3677995641, 2426754001, 1826627306, 1801558910, 1738579917,
                            4241985702, 2377822414, 2775557462, 1571967963, 306431718, 605299699,
                            3945736950, 3781665025, 269586499, 2542187119, 4141222320, 1272702028,
                            2339443181, 3839120636, 873463969, 3749228446, 2292063150, 950466071,
                            2836507061, 2230406532, 1714477543, 138451714, 3203974946, 520856374,
                            1
                        ]),
                        left: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                3808198657, 2663632093, 430236484, 4179109331, 2976821455,
                                1502387191, 2700773270, 2220243760, 3974369782, 2024520900,
                                4123080582, 2553282879, 518511851, 3723507508, 204845035, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    1833451521, 393866165, 3626261200, 3132971852, 3695593891,
                                    429046202, 2384023305, 66860
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        949559297, 1814343224, 1589915769, 16934566
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3075145729, 269665662]),
                                        left: Box::new(ProductTree::Leaf(1076187137)),
                                        right: Box::new(ProductTree::Leaf(1076211713))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[558768129, 269716982]),
                                        left: Box::new(ProductTree::Leaf(1076269057)),
                                        right: Box::new(ProductTree::Leaf(1076334593))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        816783361, 3432394685, 432106997, 16957263
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2673016833, 269797049]),
                                        left: Box::new(ProductTree::Leaf(1076391937)),
                                        right: Box::new(ProductTree::Leaf(1076531201))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[828121089, 269946949]),
                                        left: Box::new(ProductTree::Leaf(1076613121)),
                                        right: Box::new(ProductTree::Leaf(1076908033))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    1706311681, 1270083017, 4149857018, 1851761254, 1382144656,
                                    3151822842, 1934508460, 67301
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3100467201, 47214179, 1277440938, 16988510
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3378733057, 270066084]),
                                        left: Box::new(ProductTree::Leaf(1076948993)),
                                        right: Box::new(ProductTree::Leaf(1077047297))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[795475969, 270174969]),
                                        left: Box::new(ProductTree::Leaf(1077178369)),
                                        right: Box::new(ProductTree::Leaf(1077252097))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3102138369, 3032922157, 3398235019, 17014883
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2406449153, 270265380]),
                                        left: Box::new(ProductTree::Leaf(1077391361)),
                                        right: Box::new(ProductTree::Leaf(1077399553))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3648479233, 270394858]),
                                        left: Box::new(ProductTree::Leaf(1077620737)),
                                        right: Box::new(ProductTree::Leaf(1077686273))
                                    }))
                                }))
                            }))
                        })),
                        right: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                3797573633, 2277411511, 3670072177, 3710168929, 2402660037,
                                3294174062, 746428514, 2342181929, 33217290, 3563944816, 489546710,
                                3716264831, 3355120698, 624513147, 301625549, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    3994091521, 3312258343, 1969674059, 782274113, 1678305335,
                                    1980922158, 545916580, 67682
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        4076912641, 3722664394, 3096464613, 17041676
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1199374337, 270487362]),
                                        left: Box::new(ProductTree::Leaf(1077792769)),
                                        right: Box::new(ProductTree::Leaf(1077882881))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3951280129, 270598387]),
                                        left: Box::new(ProductTree::Leaf(1078038529)),
                                        right: Box::new(ProductTree::Leaf(1078079489))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3406839809, 4076482128, 3713048881, 17057741
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1736835073, 270635401]),
                                        left: Box::new(ProductTree::Leaf(1078112257)),
                                        right: Box::new(ProductTree::Leaf(1078153217))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3817488385, 270705321]),
                                        left: Box::new(ProductTree::Leaf(1078210561)),
                                        right: Box::new(ProductTree::Leaf(1078333441))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    2487836673, 1122880880, 3703095184, 4291099497, 2986154819,
                                    2754382386, 1831704171, 67914
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3642695681, 1413538466, 3554120687, 17073169
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2643320833, 270764968]),
                                        left: Box::new(ProductTree::Leaf(1078382593)),
                                        right: Box::new(ProductTree::Leaf(1078398977))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3817947137, 270820506]),
                                        left: Box::new(ProductTree::Leaf(1078497281)),
                                        right: Box::new(ProductTree::Leaf(1078505473))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3005890561, 3993276938, 4293744390, 17084714
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3079864321, 270849306]),
                                        left: Box::new(ProductTree::Leaf(1078546433)),
                                        right: Box::new(ProductTree::Leaf(1078571009))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[999768065, 270919254]),
                                        left: Box::new(ProductTree::Leaf(1078620161)),
                                        right: Box::new(ProductTree::Leaf(1078775809))
                                    }))
                                }))
                            }))
                        }))
                    }))
                })),
                right: Box::new(ProductTree::Branch(ProductBranch {
                    product: BigUint::from_slice(&[
                        1501347841, 594673369, 755763349, 263081294, 3601951008, 682262072,
                        2779995730, 623007357, 10476669, 2616168910, 30075470, 3443985448,
                        1091855876, 2805463988, 1625042599, 117805878, 3936211004, 1050396322,
                        4131165278, 3113107051, 2332973362, 2131846236, 1801625703, 2599035035,
                        3031308248, 1770964433, 2844429475, 111166169, 2985985025, 3085820881,
                        515400080, 3030949518, 2186908401, 2289970616, 2528227020, 1591941897,
                        3418105705, 3988091602, 1949119989, 3011016403, 122747760, 904571276,
                        1865685368, 2051717716, 2075887227, 238929231, 1301261873, 4216467290,
                        1208117411, 476936510, 2559465203, 2007568282, 609311345, 784524359,
                        885098663, 3923456189, 8941466, 1416259736, 3650980769, 2501958813, 1
                    ]),
                    left: Box::new(ProductTree::Branch(ProductBranch {
                        product: BigUint::from_slice(&[
                            441032705, 4165908327, 2806238415, 3341315836, 1289367134, 2204210106,
                            2378056392, 1541654266, 3863284656, 1952010702, 2902502698, 2006848139,
                            2238954672, 2833569540, 4141660701, 230029867, 3027696407, 2963073070,
                            1335563608, 2424976177, 1659813235, 3566897888, 741774183, 4224632609,
                            1450027643, 215915409, 2016170671, 278516442, 1248351561, 906059669, 1
                        ]),
                        left: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                2676219905, 1471478308, 2878197828, 4186195621, 2959469853,
                                2038460693, 871106922, 4236230574, 271914191, 3717743873,
                                3692589011, 176639771, 291417656, 2362825153, 385444153, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    446906369, 4279523643, 2976113107, 1347759102, 1122380986,
                                    3619058312, 819725731, 68286
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        4282613761, 3793812930, 1129570846, 17110940
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3952951297, 271018021]),
                                        left: Box::new(ProductTree::Leaf(1078849537)),
                                        right: Box::new(ProductTree::Leaf(1078939649))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1134968833, 271166206]),
                                        left: Box::new(ProductTree::Leaf(1079185409)),
                                        right: Box::new(ProductTree::Leaf(1079193601))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        4217356289, 1759064952, 4072471190, 17140317
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1068269569, 271269135]),
                                        left: Box::new(ProductTree::Leaf(1079357441)),
                                        right: Box::new(ProductTree::Leaf(1079431169))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[464732161, 271380321]),
                                        left: Box::new(ProductTree::Leaf(1079603201)),
                                        right: Box::new(ProductTree::Leaf(1079627777))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    2095095809, 4236349035, 2397941310, 3547416514, 3451000787,
                                    212198715, 480655727, 68541
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3312099329, 2200073215, 1591255402, 17151637
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2008317953, 271400913]),
                                        left: Box::new(ProductTree::Leaf(1079635969)),
                                        right: Box::new(ProductTree::Leaf(1079676929))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3182829569, 271427684]),
                                        left: Box::new(ProductTree::Leaf(1079685121)),
                                        right: Box::new(ProductTree::Leaf(1079734273))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        4151705601, 1834949349, 3709064829, 17163482
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[498737153, 271493589]),
                                        left: Box::new(ProductTree::Leaf(1079832577)),
                                        right: Box::new(ProductTree::Leaf(1079848961))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2646335489, 271522424]),
                                        left: Box::new(ProductTree::Leaf(1079873537)),
                                        right: Box::new(ProductTree::Leaf(1079922689))
                                    }))
                                }))
                            }))
                        })),
                        right: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                851820545, 3968090980, 3405481238, 3833604154, 3141720543,
                                2185911118, 4175927848, 2817238831, 2976528410, 336340522,
                                4241221569, 1292896281, 4238532689, 500721151, 477741463, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    1362485249, 336102286, 1030431572, 1134294127, 3271673962,
                                    2778343218, 2352612852, 68896
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2878341121, 4016239572, 2427563813, 17190580
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2579800065, 271666619]),
                                        left: Box::new(ProductTree::Leaf(1080020993)),
                                        right: Box::new(ProductTree::Leaf(1080348673))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1976262657, 271777892]),
                                        left: Box::new(ProductTree::Leaf(1080365057)),
                                        right: Box::new(ProductTree::Leaf(1080446977))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2779111425, 2037716427, 2791894081, 17213404
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[634281985, 271827352]),
                                        left: Box::new(ProductTree::Leaf(1080496129)),
                                        right: Box::new(ProductTree::Leaf(1080512513))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2950135809, 271977817]),
                                        left: Box::new(ProductTree::Leaf(1080741889)),
                                        right: Box::new(ProductTree::Leaf(1080864769))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    2643451905, 1053266254, 1955874894, 3354016831, 1056452491,
                                    1457406540, 2381744956, 69273
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1472069633, 3991762268, 4100367911, 17238602
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3185303553, 272049974]),
                                        left: Box::new(ProductTree::Leaf(1080938497)),
                                        right: Box::new(ProductTree::Leaf(1080954881))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1575100417, 272153070]),
                                        left: Box::new(ProductTree::Leaf(1081077761)),
                                        right: Box::new(ProductTree::Leaf(1081225217))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        97640449, 3450533646, 114811373, 17259383
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2078785537, 272245875]),
                                        left: Box::new(ProductTree::Leaf(1081323521)),
                                        right: Box::new(ProductTree::Leaf(1081348097))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[300556289, 272285064]),
                                        left: Box::new(ProductTree::Leaf(1081397249)),
                                        right: Box::new(ProductTree::Leaf(1081430017))
                                    }))
                                }))
                            }))
                        }))
                    })),
                    right: Box::new(ProductTree::Branch(ProductBranch {
                        product: BigUint::from_slice(&[
                            724770817, 154208670, 4112980984, 2520591591, 1112490534, 2441725477,
                            917089785, 1037005297, 4028316816, 111025289, 3082108413, 1121460348,
                            4088376772, 4193651480, 53167180, 213048063, 947144432, 2745644653,
                            2240374028, 906435690, 911070868, 2612693173, 1717648243, 1831837768,
                            4028612107, 2923484405, 273824623, 3013500514, 4253454638, 1317881002,
                            1
                        ]),
                        left: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                4126400513, 999556096, 2685871412, 2816596311, 1956616891,
                                2338760137, 3694355871, 3217436967, 3700417949, 3885351727,
                                3374256998, 3629942407, 1456079472, 2681812848, 565928432, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    63832065, 1992832695, 3899666368, 1797994890, 3283000468,
                                    3805270877, 3632707125, 69534
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1306517505, 2140879729, 2039123888, 17274031
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1038999553, 272346946]),
                                        left: Box::new(ProductTree::Leaf(1081495553)),
                                        right: Box::new(ProductTree::Leaf(1081577473))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1341259777, 272415025]),
                                        left: Box::new(ProductTree::Leaf(1081643009)),
                                        right: Box::new(ProductTree::Leaf(1081700353))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2045648897, 1019902599, 1225141167, 17288951
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2885009409, 272476922]),
                                        left: Box::new(ProductTree::Leaf(1081774081)),
                                        right: Box::new(ProductTree::Leaf(1081815041))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2784518145, 272520254]),
                                        left: Box::new(ProductTree::Leaf(1081839617)),
                                        right: Box::new(ProductTree::Leaf(1081921537))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    2988826625, 2868961583, 663792791, 388829014, 1260448571,
                                    3553708437, 3850104575, 69905
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3892740097, 3144749891, 2286360356, 17314493
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[436101121, 272619313]),
                                        left: Box::new(ProductTree::Leaf(1082060801)),
                                        right: Box::new(ProductTree::Leaf(1082093569))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[34086913, 272780320]),
                                        left: Box::new(ProductTree::Leaf(1082331137)),
                                        right: Box::new(ProductTree::Leaf(1082462209))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3659489281, 3997054239, 1371648313, 17340590
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3322781697, 272871166]),
                                        left: Box::new(ProductTree::Leaf(1082552321)),
                                        right: Box::new(ProductTree::Leaf(1082601473))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1007796225, 272939311]),
                                        left: Box::new(ProductTree::Leaf(1082699777)),
                                        right: Box::new(ProductTree::Leaf(1082724353))
                                    }))
                                }))
                            }))
                        })),
                        right: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                490684417, 4184547516, 1278199995, 34697968, 3119477143,
                                2893540450, 659145492, 308852631, 1353620567, 2624428858,
                                4068259638, 1042771712, 2977769938, 2566307446, 664406701, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    644833281, 2865706703, 3856572368, 1133061944, 1918739181,
                                    4073534986, 2553402059, 70217
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1110417409, 2424079905, 4091102741, 17357653
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1377181697, 273011594]),
                                        left: Box::new(ProductTree::Leaf(1082806273)),
                                        right: Box::new(ProductTree::Leaf(1082904577))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3222896641, 273067362]),
                                        left: Box::new(ProductTree::Leaf(1082929153)),
                                        right: Box::new(ProductTree::Leaf(1083002881))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3426729985, 239788497, 3062532863, 17374598
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[4162682881, 273133465]),
                                        left: Box::new(ProductTree::Leaf(1083076609)),
                                        right: Box::new(ProductTree::Leaf(1083117569))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1948401665, 273211973]),
                                        left: Box::new(ProductTree::Leaf(1083215873)),
                                        right: Box::new(ProductTree::Leaf(1083289601))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    3805274113, 1474675068, 3785364468, 208167774, 360772559,
                                    2999682049, 2797300493, 70628
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2354642945, 1053262983, 1396754545, 17401156
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3190300673, 273309089]),
                                        left: Box::new(ProductTree::Leaf(1083371521)),
                                        right: Box::new(ProductTree::Leaf(1083518977))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2989547521, 273453756]),
                                        left: Box::new(ProductTree::Leaf(1083535361)),
                                        right: Box::new(ProductTree::Leaf(1083928577))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1987502081, 696445729, 420692998, 17432620
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[809074689, 273596410]),
                                        left: Box::new(ProductTree::Leaf(1083953153)),
                                        right: Box::new(ProductTree::Leaf(1084076033))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[909991937, 273660510]),
                                        left: Box::new(ProductTree::Leaf(1084133377)),
                                        right: Box::new(ProductTree::Leaf(1084149761))
                                    }))
                                }))
                            }))
                        }))
                    }))
                }))
            })),
            right: Box::new(ProductTree::Branch(ProductBranch {
                product: BigUint::from_slice(&[
                    2238898177, 3714715509, 889202403, 1459927674, 3243515927, 3732391781,
                    2166457512, 605035568, 3838461149, 47616395, 2948430963, 183561900, 2170754715,
                    3535327536, 2244407623, 1329326760, 4192028561, 1993705794, 1588992594,
                    1639996008, 3192766842, 3587956349, 1389769935, 1024835679, 1394608165,
                    1629284253, 3329001819, 1658098517, 3860948953, 19387746, 673698541,
                    2784058312, 3742985405, 4076667290, 17941872, 3393682202, 3648454270, 38107654,
                    2359965693, 2222846551, 3124987175, 2234776981, 3853641484, 3510298719,
                    65308306, 3209832699, 243236028, 158800585, 2153024083, 2581555071, 3537776009,
                    1722866626, 3780464937, 533137275, 337654066, 328327431, 4016930672,
                    1587307534, 2397177796, 2388047934, 3458379569, 530782560, 2420865880,
                    2314301296, 2764023422, 3706284719, 2402561583, 3966217528, 1544485687,
                    287624492, 2033935825, 2895054198, 584756240, 4039796203, 4006869768,
                    1796531790, 4223389742, 1871763147, 932439877, 2708076574, 1668518533,
                    4068950796, 1641785674, 1575585861, 1432901990, 3121632587, 1058095064,
                    885678781, 4155600198, 2333832960, 3024860146, 573673393, 281501374,
                    3189721605, 208195436, 2690221026, 4268505041, 217241530, 3443150406,
                    3609563650, 2739128501, 14281434, 2811647256, 2077067279, 252833372,
                    2077979501, 1858470218, 1898339734, 1535949244, 3722162405, 1976060357,
                    2206004912, 17687439, 1183551236, 279845485, 1477607359, 370674952, 1883617684,
                    3382557998, 3345272972, 6
                ]),
                left: Box::new(ProductTree::Branch(ProductBranch {
                    product: BigUint::from_slice(&[
                        1974222849, 2651296343, 140112495, 3960564811, 2853473023, 3008224364,
                        3574498879, 1723405420, 3439082937, 2454015373, 3457466509, 2279437290,
                        1144909748, 991610059, 2474919399, 3739789354, 3217355952, 2900167109,
                        3866676593, 844582980, 3367028823, 3609418187, 3061275218, 4179151728,
                        2784410180, 2260890063, 504489198, 47645377, 2812642909, 689423893,
                        1698791291, 3600886012, 814449127, 3078385510, 3242519599, 4178988854,
                        1042775445, 1973416977, 3205744973, 1011204200, 3942503396, 2606080345,
                        2749927938, 1638722959, 855978164, 2329013405, 3886892069, 1208117991,
                        2199928851, 1243995366, 2513526473, 3177544500, 2026382183, 470527463,
                        1855506431, 1680101707, 787353695, 702596359, 2176216110, 699279013, 2
                    ]),
                    left: Box::new(ProductTree::Branch(ProductBranch {
                        product: BigUint::from_slice(&[
                            2287968257, 2653979533, 3969094069, 3820967002, 607565535, 3680517561,
                            1358395744, 1060061295, 644894878, 3588344962, 940010442, 1537779657,
                            2005290996, 539137302, 2919021700, 1284100011, 3278407213, 247645059,
                            3177050915, 4185441236, 2909236518, 3328551003, 1128396190, 478003683,
                            1928202693, 781565666, 1787994307, 2600715094, 3370670244, 1785786845,
                            1
                        ]),
                        left: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                916996097, 3944513136, 64389896, 2488686175, 1949425280, 171522233,
                                1332082794, 2353691456, 1906185810, 750724726, 3837950466,
                                2919234813, 2681861815, 4203765413, 773820157, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    2133647361, 3984891154, 1552206986, 2201586508, 1422386717,
                                    2645079450, 3518788624, 71026
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2056232961, 1272924161, 1903210583, 17458717
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2386878465, 273784592]),
                                        left: Box::new(ProductTree::Leaf(1084297217)),
                                        right: Box::new(ProductTree::Leaf(1084477441))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2219491329, 273881812]),
                                        left: Box::new(ProductTree::Leaf(1084518401)),
                                        right: Box::new(ProductTree::Leaf(1084641281))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        4238163969, 2805467863, 633466264, 17473097
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1548566529, 273923188]),
                                        left: Box::new(ProductTree::Leaf(1084649473)),
                                        right: Box::new(ProductTree::Leaf(1084674049))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2689597441, 273968703]),
                                        left: Box::new(ProductTree::Leaf(1084690433)),
                                        right: Box::new(ProductTree::Leaf(1084813313))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    1132158977, 4187732093, 3003098652, 2183241338, 1604534162,
                                    3329410803, 1785230579, 71364
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3971104769, 1927644069, 4200355809, 17495276
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2488647681, 274063886]),
                                        left: Box::new(ProductTree::Leaf(1084911617)),
                                        right: Box::new(ProductTree::Leaf(1084968961))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2019328001, 274175643]),
                                        left: Box::new(ProductTree::Leaf(1085140993)),
                                        right: Box::new(ProductTree::Leaf(1085181953))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3469287425, 2629316556, 2121685999, 17519461
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3697393665, 274262579]),
                                        left: Box::new(ProductTree::Leaf(1085255681)),
                                        right: Box::new(ProductTree::Leaf(1085411329))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[174546945, 274355744]),
                                        left: Box::new(ProductTree::Leaf(1085501441)),
                                        right: Box::new(ProductTree::Leaf(1085534209))
                                    }))
                                }))
                            }))
                        })),
                        right: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                2042060801, 3788986392, 2828417155, 2380722327, 3092315941,
                                1220603351, 1262784332, 1374265597, 2988403313, 2025243478,
                                957417482, 4095802783, 1426255750, 3996225598, 857476046, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    2679578625, 2652117264, 2984787435, 830906985, 3748906647,
                                    3833895471, 1413163030, 71622
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1624571905, 568134389, 1904930421, 17532026
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3026796545, 274386801]),
                                        left: Box::new(ProductTree::Leaf(1085550593)),
                                        right: Box::new(ProductTree::Leaf(1085607937))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[4234919937, 274428214]),
                                        left: Box::new(ProductTree::Leaf(1085648897)),
                                        right: Box::new(ProductTree::Leaf(1085673473))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2665619457, 1422356470, 2422620921, 17545921
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2255462401, 274492411]),
                                        left: Box::new(ProductTree::Leaf(1085779969)),
                                        right: Box::new(ProductTree::Leaf(1085796353))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[141721601, 274540046]),
                                        left: Box::new(ProductTree::Leaf(1085870081)),
                                        right: Box::new(ProductTree::Leaf(1085894657))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    302006273, 3616732890, 605934941, 971364382, 3999797957,
                                    2692126501, 276747715, 71939
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        351502337, 3456309266, 970836372, 17564328
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3866517505, 274604254]),
                                        left: Box::new(ProductTree::Leaf(1085952001)),
                                        right: Box::new(ProductTree::Leaf(1086066689))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2189238273, 274716119]),
                                        left: Box::new(ProductTree::Leaf(1086115841)),
                                        right: Box::new(ProductTree::Leaf(1086345217))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1561116673, 1794399783, 1110216666, 17591104
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3800244225, 274815578]),
                                        left: Box::new(ProductTree::Leaf(1086410753)),
                                        right: Box::new(ProductTree::Leaf(1086443521))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1720295425, 274923342]),
                                        left: Box::new(ProductTree::Leaf(1086566401)),
                                        right: Box::new(ProductTree::Leaf(1086713857))
                                    }))
                                }))
                            }))
                        }))
                    })),
                    right: Box::new(ProductTree::Branch(ProductBranch {
                        product: BigUint::from_slice(&[
                            2840371201, 3011951754, 610223932, 4129469349, 2266344050, 775179140,
                            1361035267, 2457952332, 2054870338, 3992396552, 2013546737, 1901936603,
                            4056028034, 3698122484, 1373565094, 4197248924, 2115409286, 2541577395,
                            1788298952, 4222994534, 3377012043, 2477772922, 371058088, 690957095,
                            2236752520, 709747802, 4104212716, 3881189197, 1390622516, 2266203853,
                            1
                        ]),
                        left: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                16171009, 4132472633, 3363314945, 3919557975, 1440233843,
                                543441323, 1330447791, 2359365304, 1028908963, 721623146,
                                105563469, 3820381217, 1551579380, 3577138820, 957678929, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    2253733889, 4087405598, 759034562, 1324618302, 2999546903,
                                    3008031801, 3414464664, 72308
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1461772289, 1585516972, 2642561679, 17612467
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1821294593, 275008326]),
                                        left: Box::new(ProductTree::Leaf(1086763009)),
                                        right: Box::new(ProductTree::Leaf(1086853121))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2459049985, 275064298]),
                                        left: Box::new(ProductTree::Leaf(1086902273)),
                                        right: Box::new(ProductTree::Leaf(1086935041))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1530159105, 1325123448, 2423299908, 17633185
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3298246657, 275149302]),
                                        left: Box::new(ProductTree::Leaf(1087025153)),
                                        right: Box::new(ProductTree::Leaf(1087148033))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1855791105, 275246764]),
                                        left: Box::new(ProductTree::Leaf(1087254529)),
                                        right: Box::new(ProductTree::Leaf(1087303681))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    2594275329, 4008715218, 2939585797, 2649329839, 2594380134,
                                    2376873297, 3744947124, 72641
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        4249174017, 22758294, 3492222053, 17651129
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1386340353, 275325575]),
                                        left: Box::new(ProductTree::Leaf(1087418369)),
                                        right: Box::new(ProductTree::Leaf(1087451137))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1453547521, 275350465]),
                                        left: Box::new(ProductTree::Leaf(1087475713)),
                                        right: Box::new(ProductTree::Leaf(1087492097))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1297891329, 3265890247, 3154094339, 17675608
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1890385921, 275510202]),
                                        left: Box::new(ProductTree::Leaf(1087762433)),
                                        right: Box::new(ProductTree::Leaf(1087836161))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2628730881, 275547550]),
                                        left: Box::new(ProductTree::Leaf(1087860737)),
                                        right: Box::new(ProductTree::Leaf(1087885313))
                                    }))
                                }))
                            }))
                        })),
                        right: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                5627905, 1820721804, 3251612320, 3337300550, 2884099393, 524875258,
                                918313860, 907123470, 2101505925, 3477534587, 2440942858,
                                298812490, 3126010244, 2120392301, 1069950556, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    2298994689, 2972803497, 2377551361, 4211451200, 128494584,
                                    694392753, 797543559, 73091
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1803362305, 4122914903, 3854293465, 17710641
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2864553985, 275786221]),
                                        left: Box::new(ProductTree::Leaf(1088311297)),
                                        right: Box::new(ProductTree::Leaf(1088376833))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2428469249, 275817360]),
                                        left: Box::new(ProductTree::Leaf(1088401409)),
                                        right: Box::new(ProductTree::Leaf(1088409601))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        92979201, 3701906849, 1811811090, 17725176
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[381902849, 275881718]),
                                        left: Box::new(ProductTree::Leaf(1088458753)),
                                        right: Box::new(ProductTree::Leaf(1088606209))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[180838401, 275948162]),
                                        left: Box::new(ProductTree::Leaf(1088647169)),
                                        right: Box::new(ProductTree::Leaf(1088679937))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    927858689, 3427250012, 3273915982, 1710218321, 2353981593,
                                    992885598, 1458623288, 73400
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3415777281, 2538208999, 745404832, 17739987
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[483082241, 276012536]),
                                        left: Box::new(ProductTree::Leaf(1088778241)),
                                        right: Box::new(ProductTree::Leaf(1088802817))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[919429121, 276047841]),
                                        left: Box::new(ProductTree::Leaf(1088851969)),
                                        right: Box::new(ProductTree::Leaf(1088868353))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        263544833, 2550308818, 3598833383, 17770703
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3503587329, 276166229]),
                                        left: Box::new(ProductTree::Leaf(1088966657)),
                                        right: Box::new(ProductTree::Leaf(1089220609))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1323360257, 276371922]),
                                        left: Box::new(ProductTree::Leaf(1089458177)),
                                        right: Box::new(ProductTree::Leaf(1089540097))
                                    }))
                                }))
                            }))
                        }))
                    }))
                })),
                right: Box::new(ProductTree::Branch(ProductBranch {
                    product: BigUint::from_slice(&[
                        3888553985, 441606808, 3196294177, 1718549921, 1060212768, 1829031991,
                        2479572943, 3877515798, 3126427159, 3739748061, 1910098648, 3647867357,
                        2190078556, 2667844182, 2657068709, 3407781840, 2903199403, 2501135304,
                        3484654094, 2762440626, 3259575960, 3841385922, 3291762739, 3915311017,
                        3052152128, 769272943, 3417429930, 768953331, 79009943, 544955383,
                        2128705765, 608724455, 3629709358, 2128082018, 2326001675, 4175056904,
                        2311439038, 3263194659, 1459149888, 4261881516, 3367660942, 3344479749,
                        82088527, 814839183, 636007006, 1538888827, 1779343491, 1696411119,
                        3965087343, 317809250, 1446294746, 85026981, 2381289080, 2032822381,
                        323523806, 1014371405, 1760347719, 1339418431, 2694889018, 576765349, 3
                    ]),
                    left: Box::new(ProductTree::Branch(ProductBranch {
                        product: BigUint::from_slice(&[
                            3085787137, 792512018, 2613710356, 3408217360, 2490286891, 2072324543,
                            691416708, 1157449359, 1111198970, 3669404272, 1254653691, 899589825,
                            1509857729, 1927575728, 3037699169, 2673098891, 4105648649, 645202215,
                            3012941913, 3436851484, 1670295108, 2726007917, 34203086, 3060064758,
                            2141846728, 850356247, 1878775941, 65972089, 3107393634, 2970981580, 1
                        ]),
                        left: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                3255476225, 540222920, 946538466, 2480425550, 428557182, 364350513,
                                637412234, 3615963879, 3153728206, 3793719350, 554721272,
                                2858070090, 3014148837, 1569395826, 1212738246, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    165142529, 527715100, 2911373227, 3869940248, 501904381,
                                    830942841, 2389937126, 74011
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        266493953, 115558656, 1261096605, 17818865
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[787275777, 276571451]),
                                        left: Box::new(ProductTree::Leaf(1089835009)),
                                        right: Box::new(ProductTree::Leaf(1089949697))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[284524545, 276714908]),
                                        left: Box::new(ProductTree::Leaf(1090170881)),
                                        right: Box::new(ProductTree::Leaf(1090179073))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3925180417, 3151445284, 1991034268, 17839363
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[217563137, 276752337]),
                                        left: Box::new(ProductTree::Leaf(1090203649)),
                                        right: Box::new(ProductTree::Leaf(1090293761))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1291698177, 276852161]),
                                        left: Box::new(ProductTree::Leaf(1090400257)),
                                        right: Box::new(ProductTree::Leaf(1090490369))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    271761409, 1352567073, 165226009, 3671525584, 2761045609,
                                    3662825480, 3572247679, 74416
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2114633729, 566399907, 1686865529, 17862159
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3875643393, 276916639]),
                                        left: Box::new(ProductTree::Leaf(1090498561)),
                                        right: Box::new(ProductTree::Leaf(1090646017))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3272155137, 277041460]),
                                        left: Box::new(ProductTree::Leaf(1090768897)),
                                        right: Box::new(ProductTree::Leaf(1090867201))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2452094977, 2090656335, 1326069951, 17893573
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[420511745, 277164228]),
                                        left: Box::new(ProductTree::Leaf(1091014657)),
                                        right: Box::new(ProductTree::Leaf(1091104769))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3105325057, 277280775]),
                                        left: Box::new(ProductTree::Leaf(1091178497)),
                                        right: Box::new(ProductTree::Leaf(1091399681))
                                    }))
                                }))
                            }))
                        })),
                        right: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                3185754113, 3522029888, 3180217020, 1707657375, 305106233,
                                1267529753, 3585740008, 1491884256, 340353003, 708770954,
                                1902130172, 222553437, 325658194, 3321729503, 1371096831, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    3166896129, 2700465320, 4097893512, 4044216437, 938670890,
                                    55250756, 610934374, 75064
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        4065492993, 2030933736, 3794918064, 17939290
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3341049857, 277495206]),
                                        left: Box::new(ProductTree::Leaf(1091571713)),
                                        right: Box::new(ProductTree::Leaf(1091850241))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1596858369, 277657652]),
                                        left: Box::new(ProductTree::Leaf(1092022273)),
                                        right: Box::new(ProductTree::Leaf(1092038657))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2322628609, 2647541630, 1654065990, 17971615
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1630830593, 277763889]),
                                        left: Box::new(ProductTree::Leaf(1092210689)),
                                        right: Box::new(ProductTree::Leaf(1092268033))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[691798017, 277888895]),
                                        left: Box::new(ProductTree::Leaf(1092333569)),
                                        right: Box::new(ProductTree::Leaf(1092636673))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    3240083457, 864266344, 446136622, 1386724881, 3195712175,
                                    2077016559, 4150708253, 75482
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        478617601, 333808124, 2995112351, 17996022
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3846258689, 277976424]),
                                        left: Box::new(ProductTree::Leaf(1092653057)),
                                        right: Box::new(ProductTree::Leaf(1092661249))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3745898497, 278053539]),
                                        left: Box::new(ProductTree::Leaf(1092734977)),
                                        right: Box::new(ProductTree::Leaf(1092882433))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        815308801, 645901917, 4038633656, 18014917
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3242885121, 278130667]),
                                        left: Box::new(ProductTree::Leaf(1092923393)),
                                        right: Box::new(ProductTree::Leaf(1092997121))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[525213697, 278191125]),
                                        left: Box::new(ProductTree::Leaf(1093005313)),
                                        right: Box::new(ProductTree::Leaf(1093152769))
                                    }))
                                }))
                            }))
                        }))
                    })),
                    right: Box::new(ProductTree::Branch(ProductBranch {
                        product: BigUint::from_slice(&[
                            2010726401, 2785816244, 638976528, 1211916180, 278955068, 4031081594,
                            1669327756, 1534658191, 2683158496, 3655947045, 95054458, 537955227,
                            713733654, 3645438573, 2660506785, 1539370418, 2450320536, 2138912921,
                            1420492503, 3800618973, 3316269776, 1555173629, 3023564955, 1644409320,
                            3553231259, 2446143208, 1918135238, 966960803, 2212234053, 3662344476,
                            1
                        ]),
                        left: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                1263501313, 967507923, 2966500438, 597157840, 3635186149,
                                3467986131, 2003841890, 575247459, 636057479, 245572425,
                                2706359239, 2639251630, 3472382458, 3762932623, 1496352978, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    4118380545, 1444462403, 246259441, 203122753, 2341826865,
                                    166202299, 690231355, 75890
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1185988609, 2721996097, 3118978501, 18040990
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[827752449, 278330830]),
                                        left: Box::new(ProductTree::Leaf(1093292033)),
                                        right: Box::new(ProductTree::Leaf(1093414913))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1968848897, 278393396]),
                                        left: Box::new(ProductTree::Leaf(1093439489)),
                                        right: Box::new(ProductTree::Leaf(1093513217))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        516472833, 3620961384, 2085521957, 18066954
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1264361473, 278433025]),
                                        left: Box::new(ProductTree::Leaf(1093537793)),
                                        right: Box::new(ProductTree::Leaf(1093570561))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3211534337, 278691719]),
                                        left: Box::new(ProductTree::Leaf(1093939201)),
                                        right: Box::new(ProductTree::Leaf(1094184961))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    1440088065, 3069490558, 236936291, 2191543864, 901513209,
                                    2483530242, 3747755026, 76311
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3605413889, 4040127016, 991820385, 18098911
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2171731969, 278789809]),
                                        left: Box::new(ProductTree::Leaf(1094250497)),
                                        right: Box::new(ProductTree::Leaf(1094258689))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3581165569, 278827378]),
                                        left: Box::new(ProductTree::Leaf(1094299649)),
                                        right: Box::new(ProductTree::Leaf(1094356993))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        250593281, 2640751027, 3044124178, 18109210
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1467392001, 278867038]),
                                        left: Box::new(ProductTree::Leaf(1094381569)),
                                        right: Box::new(ProductTree::Leaf(1094430721))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1132011521, 278908788]),
                                        left: Box::new(ProductTree::Leaf(1094471681)),
                                        right: Box::new(ProductTree::Leaf(1094504449))
                                    }))
                                }))
                            }))
                        })),
                        right: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                210354177, 842313434, 1012234165, 4073586591, 3431910516,
                                939875789, 2387010133, 1797474735, 1770580954, 569296764,
                                1530013358, 2345016436, 2680819028, 3584328605, 1606345738, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    3726327809, 2031106670, 1305358974, 2065966666, 3167845164,
                                    3598045193, 3924877270, 76628
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2064220161, 2458081343, 2972107459, 18137149
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[58908673, 279071641]),
                                        left: Box::new(ProductTree::Leaf(1094725633)),
                                        right: Box::new(ProductTree::Leaf(1094889473))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1803984897, 279134291]),
                                        left: Box::new(ProductTree::Leaf(1094914049)),
                                        right: Box::new(ProductTree::Leaf(1094946817))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3272720385, 3499417030, 405410584, 18146108
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[4085784577, 279159352]),
                                        left: Box::new(ProductTree::Leaf(1094963201)),
                                        right: Box::new(ProductTree::Leaf(1094995969))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2609487873, 279184415]),
                                        left: Box::new(ProductTree::Leaf(1095012353)),
                                        right: Box::new(ProductTree::Leaf(1095045121))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    242122753, 2917275688, 2008833039, 3873069840, 274180720,
                                    3709481625, 2452237530, 77011
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1462681601, 3207804864, 2755395841, 18177623
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1939136513, 279372420]),
                                        left: Box::new(ProductTree::Leaf(1095331841)),
                                        right: Box::new(ProductTree::Leaf(1095462913))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3214532609, 279455999]),
                                        left: Box::new(ProductTree::Leaf(1095536641)),
                                        right: Box::new(ProductTree::Leaf(1095585793))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1732231169, 3925416771, 2860460252, 18196117
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1167998977, 279529141]),
                                        left: Box::new(ProductTree::Leaf(1095700481)),
                                        right: Box::new(ProductTree::Leaf(1095708673))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3248586753, 279583480]),
                                        left: Box::new(ProductTree::Leaf(1095774209)),
                                        right: Box::new(ProductTree::Leaf(1095847937))
                                    }))
                                }))
                            }))
                        }))
                    }))
                }))
            }))
        })),
        right: Box::new(ProductTree::Branch(ProductBranch {
            product: BigUint::from_slice(&[
                1999552513, 1652254169, 3696155199, 3232753245, 2311799294, 3687200434, 4204907460,
                3583850546, 875907123, 3217513012, 3804788618, 3214844303, 3767323997, 528058540,
                2776898212, 1918946081, 583371664, 737805631, 1793012400, 2088144653, 4132848851,
                1492058180, 220384849, 1698343007, 4260008343, 872209833, 348190100, 1586686817,
                4264973320, 2392520834, 170721567, 165830227, 874694591, 2487565725, 2535769916,
                2391285366, 1909357812, 3293471096, 1459524156, 3015269356, 2943351032, 3259822271,
                3494972644, 3001720913, 751357014, 1390971245, 4124981213, 2027681262, 2810494009,
                2394367926, 3274581922, 4065041553, 1109543905, 865686864, 1757153366, 220168163,
                3510603972, 3868219137, 4017847243, 2027797992, 3392896857, 3091138244, 539983485,
                4275681116, 720866050, 1264847912, 2574664470, 35838487, 750409523, 4178378915,
                1667752031, 3001323739, 15013719, 3272289660, 248449288, 1215190561, 2294785827,
                542655298, 3666604529, 2026679822, 2616683011, 637001061, 7076618, 2630274388,
                410366740, 2526174557, 3885529828, 3475636145, 1929491423, 2866593811, 3503709608,
                3078201008, 2709574060, 941126038, 3506927522, 2377878575, 2586729826, 2473446302,
                711355525, 3673853682, 3383997059, 2138448226, 2959628950, 3936396339, 136478688,
                3513079839, 1270739428, 177100050, 385312817, 1147245045, 2107055278, 1776655539,
                3020921867, 2318411531, 1062255883, 821997221, 3556260217, 3029332697, 828714726,
                4228700324, 455915062, 215667736, 462283353, 1430642577, 3299184343, 2573447939,
                487172646, 2581445739, 471517768, 3691172457, 492795563, 1387589480, 3322977081,
                1544480300, 4056325000, 3367741451, 1446320585, 3853777012, 3453279516, 3184128262,
                1614709176, 3381853866, 257588837, 2689737669, 327809110, 3108207907, 468855808,
                3773072045, 3561458522, 1957195732, 433344655, 3367434678, 2318329606, 2214677497,
                3013408203, 2779605733, 3125420718, 187162830, 992744016, 1353129703, 654375346,
                1657924334, 3837558528, 4175066242, 2873207418, 1300784662, 2371787672, 1365464605,
                1160829770, 1710854429, 1922294906, 960552870, 3998551043, 2667388489, 1942830978,
                344457499, 2922828383, 545143645, 1091598245, 3004276742, 2084983022, 3037307574,
                1572302959, 307949721, 2274368107, 649732812, 894485971, 252309380, 1767868575,
                3126745457, 3951201267, 3272914444, 3661057565, 283485536, 262167389, 277884001,
                1421781892, 3010966587, 2789292972, 1676356946, 1244946531, 4266371636, 1467853042,
                3619707239, 3794151227, 192545601, 4265501499, 286629765, 2815649834, 3387333168,
                2307383196, 3378438314, 423637155, 2888053451, 2446398336, 1178161755, 1089321446,
                2109772913, 553166966, 1825282767, 1068481263, 452068242, 3504786875, 3180277811,
                2629857363, 3452335243, 3238746875, 2724908513, 4013076548, 1915925784, 3225858589,
                766132512, 1610674624, 2907217623, 1860725663, 2135876317, 121015159, 1632526370,
                3610239609, 552412994, 1916
            ]),
            left: Box::new(ProductTree::Branch(ProductBranch {
                product: BigUint::from_slice(&[
                    3949682689, 520919850, 3939403505, 1067505157, 1402843010, 3654543407,
                    352504436, 2942956639, 3136232885, 2535996796, 2524164999, 1731301746,
                    3287811921, 3459910010, 190217919, 1284119190, 1216294171, 2012696800,
                    2519759008, 2035608480, 24382808, 852484197, 1087070250, 2875703245,
                    1862188358, 178112674, 1135066865, 2224981173, 553708661, 879801483,
                    2179868810, 3745423645, 3068807011, 567081518, 2801993006, 2017383314,
                    398666017, 1616544769, 3352759713, 1700816878, 1970813791, 462398, 49969676,
                    2885849513, 884488957, 4070775174, 2665677811, 1396670183, 2628766141,
                    3438610287, 889964135, 1162781846, 3496289087, 1877745193, 1955327374,
                    320624291, 711298643, 4052561949, 3870786281, 3280394089, 2572340700,
                    2691876819, 823748912, 2506552997, 1837678554, 2594324111, 2280779798,
                    2181260224, 4034775034, 2573982567, 3637435014, 4087175431, 1354048069,
                    3448345230, 148595911, 4254012972, 3137189378, 4005846308, 4056033842,
                    478239148, 2106170773, 2089291143, 214231362, 2245538555, 155854292, 814534217,
                    2258622500, 410250449, 955425261, 4126961638, 1664942991, 2210354550,
                    2517660456, 1623031765, 1514337221, 1705406224, 475945447, 4213105818,
                    3985068358, 233565314, 1919129734, 3840504563, 4135889968, 3257151014,
                    159976887, 2257596882, 637889713, 1476574650, 4009407201, 3563701888,
                    815511478, 2387445878, 1587738703, 2005905969, 1329732416, 2157258255,
                    74890387, 1550985242, 114068616, 3430929216, 23
                ]),
                left: Box::new(ProductTree::Branch(ProductBranch {
                    product: BigUint::from_slice(&[
                        454647809, 1472260767, 816666874, 2067102670, 1225302738, 3182700009,
                        267053386, 3712785567, 3859097312, 1873522332, 2728566434, 681952989,
                        2317008205, 3711471392, 735312861, 2029115997, 69328270, 703052297,
                        975383137, 1933756472, 4732618, 1878775233, 4007383620, 3802880504,
                        663097234, 3343217553, 1562510401, 1391548228, 3394741602, 3240275392,
                        3327522244, 3379239901, 566133395, 2937399429, 518420228, 3774493838,
                        816135967, 556541707, 2487107219, 1042735688, 4172448468, 3717483072,
                        2169618552, 487825712, 2469998110, 2986349300, 1993292966, 854030172,
                        2681108133, 3236241395, 3785507135, 2629501886, 2281231927, 4198820127,
                        2084079973, 2350589521, 2431885412, 2943606802, 723663449, 917796409, 4
                    ]),
                    left: Box::new(ProductTree::Branch(ProductBranch {
                        product: BigUint::from_slice(&[
                            1078706177, 2657907712, 4104903212, 2039507726, 226880827, 22752621,
                            179531618, 715035169, 1301709611, 2806349247, 2089891886, 4274686139,
                            298980824, 614699330, 1561164990, 2680060849, 2515287443, 2377286820,
                            2653317462, 439894828, 3134465679, 1554027620, 2028936814, 3307599330,
                            1138887682, 2630766447, 2392587999, 2243208374, 1675530850, 4223058289,
                            1
                        ]),
                        left: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                3954335745, 163632596, 1082009562, 3847313139, 1263794004,
                                3926476329, 2383683511, 2484476889, 2984771151, 1999583877,
                                2409162635, 794644048, 2709550916, 2229963513, 1711963738, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    4240138241, 3694208351, 2867381647, 191803691, 3653167386,
                                    2092304959, 1791447267, 77366
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3076014081, 2514379800, 1910719378, 18222795
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2846556161, 279742350]),
                                        left: Box::new(ProductTree::Leaf(1096101889)),
                                        right: Box::new(ProductTree::Leaf(1096142849))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1303199745, 279779984]),
                                        left: Box::new(ProductTree::Leaf(1096175617)),
                                        right: Box::new(ProductTree::Leaf(1096216577))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2103648257, 745794555, 1697433394, 18234646
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3215974401, 279823892]),
                                        left: Box::new(ProductTree::Leaf(1096224769)),
                                        right: Box::new(ProductTree::Leaf(1096339457))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[632504321, 279880353]),
                                        left: Box::new(ProductTree::Leaf(1096388609)),
                                        right: Box::new(ProductTree::Leaf(1096396801))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    1324810241, 547397739, 1710520140, 2134527082, 2519224754,
                                    2320743748, 2646811892, 77642
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        4218585089, 3712785937, 417699995, 18251411
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2444689409, 279943092]),
                                        left: Box::new(ProductTree::Leaf(1096486913)),
                                        right: Box::new(ProductTree::Leaf(1096544257))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3921379329, 280018388]),
                                        left: Box::new(ProductTree::Leaf(1096609793)),
                                        right: Box::new(ProductTree::Leaf(1096716289))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        327450625, 914677270, 3853789482, 18271052
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3351273473, 280099971]),
                                        left: Box::new(ProductTree::Leaf(1096765441)),
                                        right: Box::new(ProductTree::Leaf(1096880129))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3217301505, 280162736]),
                                        left: Box::new(ProductTree::Leaf(1096937473)),
                                        right: Box::new(ProductTree::Leaf(1096953857))
                                    }))
                                }))
                            }))
                        })),
                        right: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                278487041, 3261391405, 4277403678, 2937425538, 1682221070,
                                1749117882, 3334853572, 1256768445, 3042508948, 2601210893,
                                2252181011, 3689778353, 551582239, 2986344953, 1795437455, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    254656513, 1739370786, 730922512, 136652770, 1452056228,
                                    463056848, 2567419, 77897
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3314655233, 1296775301, 3933881346, 18285384
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3888537601, 280200397]),
                                        left: Box::new(ProductTree::Leaf(1096962049)),
                                        right: Box::new(ProductTree::Leaf(1097076737))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[634077185, 280282008]),
                                        left: Box::new(ProductTree::Leaf(1097175041)),
                                        right: Box::new(ProductTree::Leaf(1097183233))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2510036993, 173742348, 3596709140, 18296856
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[4123836417, 280307120]),
                                        left: Box::new(ProductTree::Leaf(1097199617)),
                                        right: Box::new(ProductTree::Leaf(1097256961))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2546950145, 280351071]),
                                        left: Box::new(ProductTree::Leaf(1097306113)),
                                        right: Box::new(ProductTree::Leaf(1097322497))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    694919169, 4073847727, 3852461233, 3677982355, 1927774141,
                                    4077123486, 1530586702, 78185
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        799817729, 1558549795, 2034295525, 18314484
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1070850049, 280426423]),
                                        left: Box::new(ProductTree::Leaf(1097453569)),
                                        right: Box::new(ProductTree::Leaf(1097469953))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2077777921, 280501783]),
                                        left: Box::new(ProductTree::Leaf(1097527297)),
                                        right: Box::new(ProductTree::Leaf(1097691137))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2713673729, 1205178087, 818562462, 18335408
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[870064129, 280564594]),
                                        left: Box::new(ProductTree::Leaf(1097715713)),
                                        right: Box::new(ProductTree::Leaf(1097748481))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[98779137, 280683950]),
                                        left: Box::new(ProductTree::Leaf(1097895937)),
                                        right: Box::new(ProductTree::Leaf(1098035201))
                                    }))
                                }))
                            }))
                        }))
                    })),
                    right: Box::new(ProductTree::Branch(ProductBranch {
                        product: BigUint::from_slice(&[
                            3268255745, 162589254, 1668103774, 1683673103, 2032257947, 2101045980,
                            2512747276, 2062820502, 921561966, 4013650488, 3918152990, 685970543,
                            1863190431, 2099001300, 1046618066, 1864232662, 4249966397, 2488884226,
                            4099066038, 4011171567, 518352103, 4151350392, 1984614094, 2482167385,
                            3574505955, 3803066059, 3983483251, 914202629, 1568214074, 535288275,
                            2
                        ]),
                        left: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                432103425, 3691610546, 4125388729, 3442101738, 2050927339,
                                4057881725, 451238256, 1962177161, 332926391, 430967610, 820396891,
                                2046089451, 3774542997, 3789248352, 1903988972, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    1404944385, 4261730436, 4048187659, 881486964, 2668069162,
                                    1672421491, 2298501605, 78569
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1205239809, 3206376387, 2485180625, 18360732
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[266936321, 280782386]),
                                        left: Box::new(ProductTree::Leaf(1098084353)),
                                        right: Box::new(ProductTree::Leaf(1098231809))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[736976897, 280853607]),
                                        left: Box::new(ProductTree::Leaf(1098280961)),
                                        right: Box::new(ProductTree::Leaf(1098313729))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3957800961, 3933364413, 2126764632, 18379091
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2280677377, 280903885]),
                                        left: Box::new(ProductTree::Leaf(1098387457)),
                                        right: Box::new(ProductTree::Leaf(1098403841))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3489062913, 281012833]),
                                        left: Box::new(ProductTree::Leaf(1098485761)),
                                        right: Box::new(ProductTree::Leaf(1098731521))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    3993214977, 2701906205, 2360018510, 3298667401, 1638274702,
                                    248226453, 3042258449, 78897
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3724107777, 1747411168, 3863677971, 18398973
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[4261134337, 281094564]),
                                        left: Box::new(ProductTree::Leaf(1098756097)),
                                        right: Box::new(ProductTree::Leaf(1098780673))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3355287553, 281126001]),
                                        left: Box::new(ProductTree::Leaf(1098805249)),
                                        right: Box::new(ProductTree::Leaf(1098854401))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        4161421313, 3420059648, 2804337217, 18417498
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3858907137, 281203552]),
                                        left: Box::new(ProductTree::Leaf(1098919937)),
                                        right: Box::new(ProductTree::Leaf(1099042817))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[302514177, 281299982]),
                                        left: Box::new(ProductTree::Leaf(1099067393)),
                                        right: Box::new(ProductTree::Leaf(1099272193))
                                    }))
                                }))
                            }))
                        })),
                        right: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                3507240961, 566537651, 3218969584, 2682500806, 2393762386,
                                720799655, 2182025591, 4029287354, 501344921, 2460356652,
                                257570354, 4014943028, 1037070787, 142654931, 2027473464, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    3528810497, 2147566232, 4247213374, 3432215489, 1223014317,
                                    1911048317, 722475438, 79283
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3525255169, 2567408214, 2118409509, 18440434
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1913454593, 281383850]),
                                        left: Box::new(ProductTree::Leaf(1099296769)),
                                        right: Box::new(ProductTree::Leaf(1099370497))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3222413313, 281469824]),
                                        left: Box::new(ProductTree::Leaf(1099411457)),
                                        right: Box::new(ProductTree::Leaf(1099591681))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2016821249, 391707177, 892045754, 18465867
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[437813249, 281576790]),
                                        left: Box::new(ProductTree::Leaf(1099665409)),
                                        right: Box::new(ProductTree::Leaf(1099755521))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2115878913, 281664890]),
                                        left: Box::new(ProductTree::Leaf(1099763713)),
                                        right: Box::new(ProductTree::Leaf(1100001281))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    179757057, 1785697438, 1627988097, 16070103, 1933912700,
                                    3503902347, 244359906, 79745
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2823995393, 589377953, 2696335824, 18497247
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1780875265, 281803369]),
                                        left: Box::new(ProductTree::Leaf(1100132353)),
                                        right: Box::new(ProductTree::Leaf(1100173313))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3190603777, 281916691]),
                                        left: Box::new(ProductTree::Leaf(1100296193)),
                                        right: Box::new(ProductTree::Leaf(1100451841))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        174333953, 2447110687, 2878919721, 18516398
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2452652033, 281979660]),
                                        left: Box::new(ProductTree::Leaf(1100492801)),
                                        right: Box::new(ProductTree::Leaf(1100500993))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2821955585, 282032138]),
                                        left: Box::new(ProductTree::Leaf(1100574721)),
                                        right: Box::new(ProductTree::Leaf(1100623873))
                                    }))
                                }))
                            }))
                        }))
                    }))
                })),
                right: Box::new(ProductTree::Branch(ProductBranch {
                    product: BigUint::from_slice(&[
                        273809409, 3532147819, 4294623418, 2519376425, 2697331265, 880028310,
                        2277867186, 1840033020, 804816049, 3233543922, 2564556452, 1458470330,
                        2664247069, 899470204, 869170569, 892887792, 2200095390, 3277364104,
                        755186688, 365239329, 3425645088, 22812069, 3735686930, 1116637044,
                        78525696, 155690140, 1611223266, 3481283063, 2221689103, 3051002980,
                        1736005893, 915751445, 3749395739, 2821677089, 442049613, 100759496,
                        1855805801, 189105402, 2646936917, 1465792662, 2773785256, 3406260828,
                        2005479421, 1186482798, 134598803, 3563768105, 1932227373, 3176984317,
                        3686632903, 2551401032, 1800892971, 2117797567, 3948118118, 1466383794,
                        1933100996, 3632668897, 2078622891, 505973219, 2203113985, 2783034802, 5
                    ]),
                    left: Box::new(ProductTree::Branch(ProductBranch {
                        product: BigUint::from_slice(&[
                            3712253953, 1579254075, 3388004944, 256898584, 3025745714, 3257173807,
                            21713218, 349042226, 676432746, 2450454759, 1135629563, 1452330251,
                            2048095843, 4251008187, 100603548, 334515762, 195112611, 295838945,
                            4061997431, 3696762525, 2325440453, 3499205889, 662303854, 2531135298,
                            3582602844, 78306455, 1101638157, 352247266, 3078476540, 1177993581, 2
                        ]),
                        left: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                3961413633, 858045916, 116501115, 2232085410, 2553680867,
                                2668877956, 3614114829, 1359983965, 1025291214, 810921430,
                                738679706, 3810581604, 249383503, 901453722, 2131529310, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    753958913, 3676322820, 1923133989, 1794146589, 2873433503,
                                    24802242, 638378196, 80019
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2020655105, 1261660990, 844480113, 18530323
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[305561601, 282080423]),
                                        left: Box::new(ProductTree::Leaf(1100689409)),
                                        right: Box::new(ProductTree::Leaf(1100697601))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3325706241, 282143408]),
                                        left: Box::new(ProductTree::Leaf(1100812289)),
                                        right: Box::new(ProductTree::Leaf(1100820481))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3833577473, 1596656513, 310318448, 18546877
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[876355585, 282174904]),
                                        left: Box::new(ProductTree::Leaf(1100861441)),
                                        right: Box::new(ProductTree::Leaf(1100894209))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3158548481, 282300903]),
                                        left: Box::new(ProductTree::Leaf(1101107201)),
                                        right: Box::new(ProductTree::Leaf(1101139969))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    3475890177, 4043242687, 1649135981, 3016166861, 3434272318,
                                    1592784960, 4227397626, 80311
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1620099073, 2410933071, 948941810, 18565652
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2051407873, 282340809]),
                                        left: Box::new(ProductTree::Leaf(1101189121)),
                                        right: Box::new(ProductTree::Leaf(1101213697))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1716174849, 282420629]),
                                        left: Box::new(ProductTree::Leaf(1101352961)),
                                        right: Box::new(ProductTree::Leaf(1101361153))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2526879745, 789322776, 869803587, 18579328
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2823577601, 282447938]),
                                        left: Box::new(ProductTree::Leaf(1101385729)),
                                        right: Box::new(ProductTree::Leaf(1101434881))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2790309889, 282521470]),
                                        left: Box::new(ProductTree::Leaf(1101549569)),
                                        right: Box::new(ProductTree::Leaf(1101557761))
                                    }))
                                }))
                            }))
                        })),
                        right: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                3844481025, 2552360312, 2411451668, 1187585056, 2342914792,
                                2591611771, 11702454, 3068169279, 3295205642, 3743046612,
                                2840261657, 2887837105, 2437939512, 1502980502, 2233151308, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    1701068801, 1704358564, 3920840579, 2414745379, 2446799261,
                                    2544657115, 843778006, 80572
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1621606401, 1301696600, 1313104442, 18591076
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[173178881, 282550886]),
                                        left: Box::new(ProductTree::Leaf(1101598721)),
                                        right: Box::new(ProductTree::Leaf(1101623297))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[911556609, 282597113]),
                                        left: Box::new(ProductTree::Leaf(1101672449)),
                                        right: Box::new(ProductTree::Leaf(1101729793))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        616333313, 3467583431, 856811834, 18614035
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1918566401, 282693782]),
                                        left: Box::new(ProductTree::Leaf(1101852673)),
                                        right: Box::new(ProductTree::Leaf(1101926401))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[845250561, 282803080]),
                                        left: Box::new(ProductTree::Leaf(1102041089)),
                                        right: Box::new(ProductTree::Leaf(1102163969))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    2277629953, 2009158600, 3086909481, 2292976916, 2910174723,
                                    1113231191, 4193845103, 81021
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2161524737, 2920592154, 2160465194, 18642556
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3429392385, 282918707]),
                                        left: Box::new(ProductTree::Leaf(1102295041)),
                                        right: Box::new(ProductTree::Leaf(1102360577))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2154684417, 283011226]),
                                        left: Box::new(ProductTree::Leaf(1102483457)),
                                        right: Box::new(ProductTree::Leaf(1102532609))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1592500225, 3779616763, 1602198559, 18666256
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[544333825, 283078521]),
                                        left: Box::new(ProductTree::Leaf(1102565377)),
                                        right: Box::new(ProductTree::Leaf(1102712833))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[4135174145, 283211033]),
                                        left: Box::new(ProductTree::Leaf(1102860289)),
                                        right: Box::new(ProductTree::Leaf(1102934017))
                                    }))
                                }))
                            }))
                        }))
                    })),
                    right: Box::new(ProductTree::Branch(ProductBranch {
                        product: BigUint::from_slice(&[
                            1930264577, 961010322, 3398997384, 3863470969, 2220576623, 899490284,
                            2301634153, 2085576016, 1829946483, 1445657827, 2543337027, 2712578916,
                            1701646563, 543104939, 605702863, 466385217, 689437311, 3902591001,
                            3255665598, 1864076274, 989591162, 3186334566, 1387360105, 2449467936,
                            1894239167, 3647589282, 813925410, 1182399517, 585048346, 2076274452,
                            2
                        ]),
                        left: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                1959755777, 941636148, 3266056023, 3052930544, 713366839,
                                2310642667, 374644655, 719121523, 1222426504, 959317230, 533830925,
                                128737589, 2190551728, 1742639534, 2400230312, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    36814849, 3895159194, 915122944, 4287241491, 3473279862,
                                    2825793390, 1293696296, 81562
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2433220609, 1319074966, 1761310261, 18697753
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2156011521, 283351993]),
                                        left: Box::new(ProductTree::Leaf(1103147009)),
                                        right: Box::new(ProductTree::Leaf(1103196161))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2491801601, 283415120]),
                                        left: Box::new(ProductTree::Leaf(1103245313)),
                                        right: Box::new(ProductTree::Leaf(1103343617))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1898561537, 3713089013, 3182449690, 18735267
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3130130433, 283619276]),
                                        left: Box::new(ProductTree::Leaf(1103548417)),
                                        right: Box::new(ProductTree::Leaf(1103835137))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[915914753, 283716125]),
                                        left: Box::new(ProductTree::Leaf(1103843329)),
                                        right: Box::new(ProductTree::Leaf(1103917057))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    1654505473, 468675314, 2138558712, 2059392723, 472999837,
                                    2113579962, 3925766646, 82086
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3007356929, 469242807, 382837970, 18760726
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2862284801, 283770870]),
                                        left: Box::new(ProductTree::Leaf(1103966209)),
                                        right: Box::new(ProductTree::Leaf(1104007169))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[681943041, 283949881]),
                                        left: Box::new(ProductTree::Leaf(1104261121)),
                                        right: Box::new(ProductTree::Leaf(1104408577))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2069700609, 2302227832, 890341586, 18792482
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2360107009, 284063634]),
                                        left: Box::new(ProductTree::Leaf(1104457729)),
                                        right: Box::new(ProductTree::Leaf(1104654337))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2796601345, 284137378]),
                                        left: Box::new(ProductTree::Leaf(1104695297)),
                                        right: Box::new(ProductTree::Leaf(1104703489))
                                    }))
                                }))
                            }))
                        })),
                        right: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                4265476097, 3227219228, 1352366016, 2041003025, 3755938395,
                                1221335109, 3498986500, 1937782558, 4076141078, 1428069587,
                                434563239, 2683981555, 3706172501, 3053283582, 2547402667, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    3506659329, 1014851312, 3089094883, 3011137277, 1705993372,
                                    1021920290, 1082550682, 82583
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3581231105, 2070397401, 3947882683, 18819392
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2293473281, 284185841]),
                                        left: Box::new(ProductTree::Leaf(1104719873)),
                                        right: Box::new(ProductTree::Leaf(1104867329))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1757519873, 284421900]),
                                        left: Box::new(ProductTree::Leaf(1105195009)),
                                        right: Box::new(ProductTree::Leaf(1105309697))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3549306881, 3589857182, 2887748397, 18847173
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2932195329, 284491472]),
                                        left: Box::new(ProductTree::Leaf(1105367041)),
                                        right: Box::new(ProductTree::Leaf(1105408001))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[4240990209, 284535749]),
                                        left: Box::new(ProductTree::Leaf(1105457153)),
                                        right: Box::new(ProductTree::Leaf(1105489921))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    1362796545, 2262687590, 4049451236, 1104459386, 1030056579,
                                    1428949499, 894443447, 82854
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        697712641, 53308834, 3990641500, 18856253
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2932490241, 284567378]),
                                        left: Box::new(ProductTree::Leaf(1105514497)),
                                        right: Box::new(ProductTree::Leaf(1105555457))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2060189697, 284596900]),
                                        left: Box::new(ProductTree::Leaf(1105580033)),
                                        right: Box::new(ProductTree::Leaf(1105604609))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3483656193, 3429498922, 1000269765, 18872047
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[4241514497, 284670710]),
                                        left: Box::new(ProductTree::Leaf(1105686529)),
                                        right: Box::new(ProductTree::Leaf(1105784833))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1389625345, 284731876]),
                                        left: Box::new(ProductTree::Leaf(1105825793)),
                                        right: Box::new(ProductTree::Leaf(1105883137))
                                    }))
                                }))
                            }))
                        }))
                    }))
                }))
            })),
            right: Box::new(ProductTree::Branch(ProductBranch {
                product: BigUint::from_slice(&[
                    667115521, 1037927420, 950538018, 773383260, 3689527709, 1680099691,
                    3780292881, 1634705669, 820871321, 1468452201, 1500443759, 2002971869,
                    3023611186, 3210245112, 1841996860, 1012724082, 897044505, 4016681255,
                    649112032, 3121286964, 2924394367, 1459516640, 1971066647, 4054870497,
                    3221936019, 3780990709, 1530503090, 4015341299, 2563126001, 1811453741,
                    865285032, 1033524536, 4268950505, 2991781899, 2821948078, 3333626706,
                    2720092072, 2893759960, 2437471737, 1137382275, 2459866071, 2908326879,
                    3346243000, 4004483707, 4093651340, 2668996764, 829011054, 1140208563,
                    1045727328, 1721233165, 365479796, 2168707153, 889251067, 1856492601,
                    4224756705, 3764571502, 4148073953, 3044690375, 89132313, 3439019770,
                    2782321982, 1574850962, 1160334229, 3281028968, 928602130, 3533866118,
                    2987163992, 3135732310, 68075023, 3434773951, 2630701447, 310131894,
                    4127284740, 3887440962, 2102020072, 908499043, 1106511774, 349968702,
                    1075358774, 2534961526, 2595826644, 2257146624, 1122640726, 2309029336,
                    3148859756, 459642210, 1352528149, 3848344654, 3128904129, 822439153,
                    2833107412, 1523812014, 2472982478, 626810364, 736849337, 1038515412,
                    1482602153, 3629636801, 3217522733, 2001349322, 4027660104, 2142812193,
                    360608411, 3500882080, 493688197, 4267650146, 2616853518, 3382476644,
                    2547648415, 2519787297, 730087959, 322475005, 1148487413, 419289878,
                    1466967883, 4129183182, 186908941, 1213040671, 1097989984, 2205805925, 80
                ]),
                left: Box::new(ProductTree::Branch(ProductBranch {
                    product: BigUint::from_slice(&[
                        3031433217, 4006511609, 206829889, 399981524, 1060438715, 3787163773,
                        2725874087, 1479023843, 1795380651, 3146215110, 2013147967, 949666229,
                        2839063144, 973416231, 1764614782, 304302634, 1371156112, 149266816,
                        3612528401, 771494368, 315640366, 1993802750, 312851425, 181498236,
                        695275014, 662047114, 374845910, 3736536009, 1175305669, 3698127973,
                        1857366592, 2170234694, 1236946247, 3332694581, 3276973446, 767263182,
                        2928097899, 4291284024, 2520132557, 143817733, 205244808, 832518857,
                        2555258607, 1004645787, 835083403, 1992386129, 409504430, 266279861,
                        3328812899, 558791990, 379389826, 1001986538, 2180489675, 1964240773,
                        633821104, 3144644743, 481373063, 874977973, 2958361273, 2976205428, 7
                    ]),
                    left: Box::new(ProductTree::Branch(ProductBranch {
                        product: BigUint::from_slice(&[
                            131596289, 4052483983, 1931552644, 2310515528, 100986818, 868490289,
                            3866643691, 622803722, 2336595454, 2252275558, 1522744839, 3416433769,
                            134100668, 4130142377, 3791112526, 3248542973, 1364785613, 1261229298,
                            746052750, 261455181, 2702562279, 3978881925, 2821074493, 261520944,
                            1935458268, 2528993997, 1124476182, 1278589764, 1093438565, 2880344014,
                            2
                        ]),
                        left: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                1696915457, 1894121820, 526574676, 2506957255, 399363602,
                                3170837037, 2999472053, 1684365666, 1912366696, 1491839424,
                                3463704576, 405846877, 4003258392, 3846644297, 2641659200, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    259309569, 1322227019, 1073967786, 4101349978, 4046949498,
                                    3419545684, 362920647, 83140
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2175844353, 933737202, 639374275, 18885892
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2530656257, 284778281]),
                                        left: Box::new(ProductTree::Leaf(1105907713)),
                                        right: Box::new(ProductTree::Leaf(1105981441))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1792671745, 284833129]),
                                        left: Box::new(ProductTree::Leaf(1106030593)),
                                        right: Box::new(ProductTree::Leaf(1106071553))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2378432513, 3249175877, 3320569890, 18907443
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1122009089, 284942840]),
                                        left: Box::new(ProductTree::Leaf(1106227201)),
                                        right: Box::new(ProductTree::Leaf(1106300929))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2330165249, 284993483]),
                                        left: Box::new(ProductTree::Leaf(1106350081)),
                                        right: Box::new(ProductTree::Leaf(1106374657))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    3182436353, 3887773726, 1081602533, 4139739467, 2825472114,
                                    3520355642, 4286024357, 83432
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2278432769, 2942760004, 3301563275, 18918786
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[753246209, 285029358]),
                                        left: Box::new(ProductTree::Leaf(1106399233)),
                                        right: Box::new(ProductTree::Leaf(1106464769))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2330492929, 285077898]),
                                        left: Box::new(ProductTree::Leaf(1106513921)),
                                        right: Box::new(ProductTree::Leaf(1106538497))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        4125229057, 3643112257, 1961068831, 18941066
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3169558529, 285130662]),
                                        left: Box::new(ProductTree::Leaf(1106546689)),
                                        right: Box::new(ProductTree::Leaf(1106710529))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[485908481, 285312215]),
                                        left: Box::new(ProductTree::Leaf(1106931713)),
                                        right: Box::new(ProductTree::Leaf(1107030017))
                                    }))
                                }))
                            }))
                        })),
                        right: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                1723015169, 2722801541, 2544480849, 2812884368, 1470095104,
                                1860582772, 3971261498, 1528693194, 4199759421, 1676385134,
                                3862061202, 2910990462, 2229563420, 4143209429, 2807112008, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    808714241, 3144520900, 2350838679, 1582131004, 2000061191,
                                    3797010821, 3444067447, 84086
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3356123137, 3261158125, 1566353408, 18986388
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3607085057, 285470579]),
                                        left: Box::new(ProductTree::Leaf(1107054593)),
                                        right: Box::new(ProductTree::Leaf(1107521537))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2567610369, 285654364]),
                                        left: Box::new(ProductTree::Leaf(1107619841)),
                                        right: Box::new(ProductTree::Leaf(1107668993))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2351538177, 1338158758, 719390598, 19021525
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1427300353, 285793815]),
                                        left: Box::new(ProductTree::Leaf(1107816449)),
                                        right: Box::new(ProductTree::Leaf(1108013057))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3474374657, 285859329]),
                                        left: Box::new(ProductTree::Leaf(1108021249)),
                                        right: Box::new(ProductTree::Leaf(1108062209))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    3061784577, 2488511214, 1783471336, 245889384, 665885287,
                                    4225740776, 1222342538, 84461
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1278763009, 386957956, 3922723752, 19038125
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1192869889, 285910054]),
                                        left: Box::new(ProductTree::Leaf(1108135937)),
                                        right: Box::new(ProductTree::Leaf(1108144129))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1830723585, 285992489]),
                                        left: Box::new(ProductTree::Leaf(1108217857)),
                                        right: Box::new(ProductTree::Leaf(1108381697))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3796287489, 2888397589, 1571909016, 19054315
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1730330625, 286062254]),
                                        left: Box::new(ProductTree::Leaf(1108430849)),
                                        right: Box::new(ProductTree::Leaf(1108439041))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2065956865, 286083396]),
                                        left: Box::new(ProductTree::Leaf(1108463617)),
                                        right: Box::new(ProductTree::Leaf(1108488193))
                                    }))
                                }))
                            }))
                        }))
                    })),
                    right: Box::new(ProductTree::Branch(ProductBranch {
                        product: BigUint::from_slice(&[
                            2899836929, 3923691050, 3832545552, 479829481, 21059347, 2871667169,
                            3588956400, 3056109279, 777527057, 654262232, 583838640, 2983653560,
                            2837101686, 2901994895, 3485664077, 3418764908, 1141684226, 726154504,
                            1874331750, 996912364, 3312848407, 2566292949, 3554280554, 3646980402,
                            2922983440, 4208194570, 237563956, 4267012249, 3059751912, 3782032852,
                            2
                        ]),
                        left: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                1974681601, 207829672, 1321373904, 2926458712, 3143944515,
                                201068949, 1256082, 3646342442, 3888482422, 2996371852, 375020027,
                                1230065714, 427968257, 1061983189, 2915468774, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    3904126977, 1272813101, 9814470, 1898130668, 3773585132,
                                    3164408923, 4138268655, 84726
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1683038209, 2115670850, 4190136139, 19066005
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1160175617, 286132026]),
                                        left: Box::new(ProductTree::Leaf(1108553729)),
                                        right: Box::new(ProductTree::Leaf(1108586497))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2804563969, 286189117]),
                                        left: Box::new(ProductTree::Leaf(1108611073)),
                                        right: Box::new(ProductTree::Leaf(1108750337))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        6496257, 799779586, 1908383705, 19086301
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1496276993, 286275824]),
                                        left: Box::new(ProductTree::Leaf(1108774913)),
                                        right: Box::new(ProductTree::Leaf(1108922369))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[4281581569, 286349853]),
                                        left: Box::new(ProductTree::Leaf(1108979713)),
                                        right: Box::new(ProductTree::Leaf(1109004289))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    3439263745, 176666094, 3448345569, 3191332976, 2013647493,
                                    4278181859, 102454811, 85102
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        41107457, 3551457203, 2872956743, 19104496
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[221732865, 286411199]),
                                        left: Box::new(ProductTree::Leaf(1109094401)),
                                        right: Box::new(ProductTree::Leaf(1109127169))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1631313921, 286487360]),
                                        left: Box::new(ProductTree::Leaf(1109217281)),
                                        right: Box::new(ProductTree::Leaf(1109299201))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2592849921, 3744439121, 1722736543, 19132166
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3712114689, 286597384]),
                                        left: Box::new(ProductTree::Leaf(1109323777)),
                                        right: Box::new(ProductTree::Leaf(1109618689))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2034851841, 286715906]),
                                        left: Box::new(ProductTree::Leaf(1109692417)),
                                        right: Box::new(ProductTree::Leaf(1109708801))
                                    }))
                                }))
                            }))
                        })),
                        right: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                3341074433, 2659897158, 3347263756, 2746370833, 333558994,
                                3834285984, 446365500, 2997384913, 2885268802, 1055731888,
                                2472144296, 1737769755, 4060197866, 1290703060, 3074517024, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    1163624449, 1462628045, 2839124454, 1472824361, 4168236953,
                                    2427539153, 750560738, 85568
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2326151169, 3288671160, 3055874400, 19162129
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[558899201, 286830213]),
                                        left: Box::new(ProductTree::Leaf(1109880833)),
                                        right: Box::new(ProductTree::Leaf(1109962753))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2773884929, 286931838]),
                                        left: Box::new(ProductTree::Leaf(1110085633)),
                                        right: Box::new(ProductTree::Leaf(1110151169))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3937746945, 1876093535, 3224087267, 19179105
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1062846465, 286993246]),
                                        left: Box::new(ProductTree::Leaf(1110224897)),
                                        right: Box::new(ProductTree::Leaf(1110249473))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2338029569, 287022893]),
                                        left: Box::new(ProductTree::Leaf(1110282241)),
                                        right: Box::new(ProductTree::Leaf(1110306817))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    566837249, 905891531, 3747206563, 3193996031, 1680839579,
                                    1202782913, 543988593, 86124
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2094637057, 3026507097, 405420771, 19220316
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1298423809, 287173267]),
                                        left: Box::new(ProductTree::Leaf(1110454273)),
                                        right: Box::new(ProductTree::Leaf(1110716417))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[930430977, 287459309]),
                                        left: Box::new(ProductTree::Leaf(1111109633)),
                                        right: Box::new(ProductTree::Leaf(1111166977))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1089445889, 1109718077, 27959130, 19245277
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[628563969, 287491100]),
                                        left: Box::new(ProductTree::Leaf(1111183361)),
                                        right: Box::new(ProductTree::Leaf(1111216129))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2339930113, 287514414]),
                                        left: Box::new(ProductTree::Leaf(1111232513)),
                                        right: Box::new(ProductTree::Leaf(1111257089))
                                    }))
                                }))
                            }))
                        }))
                    }))
                })),
                right: Box::new(ProductTree::Branch(ProductBranch {
                    product: BigUint::from_slice(&[
                        1930649601, 4108182320, 4019713979, 3259307767, 1390343128, 3194533210,
                        322104644, 2699570787, 3126256384, 2936201135, 3754138755, 577932894,
                        118215682, 670198465, 3673666552, 1712773965, 1571119549, 4221264669,
                        1667678606, 3189142066, 1231056698, 3939107520, 1113446618, 1931577932,
                        2318038531, 480763781, 1060157943, 410007181, 3075152280, 506660984,
                        1234015596, 4098932753, 2478041716, 3586750072, 2909982185, 1853340301,
                        1855066841, 1577874781, 2980553992, 3788823692, 3282134984, 2761248438,
                        3024440727, 2268922970, 192584377, 2606790904, 503118937, 3063478554,
                        4260915774, 1799727271, 1091445823, 2935006316, 4227089758, 1958531687,
                        3713521629, 1601316792, 4027584829, 1156026962, 1136935641, 2000977641, 10
                    ]),
                    left: Box::new(ProductTree::Branch(ProductBranch {
                        product: BigUint::from_slice(&[
                            463904769, 3246478260, 1596726875, 613039075, 1224273856, 974157667,
                            2709684262, 999778796, 1502346843, 4136603688, 3490394136, 4084797246,
                            2539831975, 3980464094, 1708478557, 1734960476, 1286897139, 2127173817,
                            813713196, 2799089885, 2947670584, 1623935954, 3464674526, 3760934406,
                            1342645661, 2488514463, 3182478885, 3837435504, 2744798100, 418545326,
                            3
                        ]),
                        left: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                910417921, 3779802279, 1662805506, 978081811, 2328495503,
                                1530061463, 940533441, 469583708, 3071162711, 1599027362,
                                2270091070, 2992705988, 698725135, 2330919218, 3200385625, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    537436161, 328098579, 4179128831, 4058792680, 1760847467,
                                    476680473, 2686245169, 86446
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3103760385, 3111479511, 2876848857, 19263443
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[159277057, 287614041]),
                                        left: Box::new(ProductTree::Leaf(1111412737)),
                                        right: Box::new(ProductTree::Leaf(1111461889))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3078701057, 287662800]),
                                        left: Box::new(ProductTree::Leaf(1111502849)),
                                        right: Box::new(ProductTree::Leaf(1111560193))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3339255809, 2368209886, 1128040502, 19274094
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2072215553, 287700963]),
                                        left: Box::new(ProductTree::Leaf(1111584769)),
                                        right: Box::new(ProductTree::Leaf(1111625729))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[193298433, 287734888]),
                                        left: Box::new(ProductTree::Leaf(1111658497)),
                                        right: Box::new(ProductTree::Leaf(1111683073))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    238764033, 4091678421, 4068831676, 2539487503, 1030435824,
                                    1378428079, 4208654752, 86704
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3205791745, 1846691964, 3765246343, 19287164
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2408103937, 287790019]),
                                        left: Box::new(ProductTree::Leaf(1111756801)),
                                        right: Box::new(ProductTree::Leaf(1111797761))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1602994177, 287840914]),
                                        left: Box::new(ProductTree::Leaf(1111822337)),
                                        right: Box::new(ProductTree::Leaf(1111928833))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        4012294145, 1988572797, 3701403965, 19307920
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1636835329, 287915145]),
                                        left: Box::new(ProductTree::Leaf(1111994369)),
                                        right: Box::new(ProductTree::Leaf(1112043521))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3113656321, 288025447]),
                                        left: Box::new(ProductTree::Leaf(1112174593)),
                                        right: Box::new(ProductTree::Leaf(1112289281))
                                    }))
                                }))
                            }))
                        })),
                        right: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                1029881857, 3318711294, 2824133147, 91569589, 254598887,
                                3749515273, 2691489048, 398068225, 957189363, 4098032317,
                                1705968586, 522047603, 2528153554, 3238156692, 3328138821, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    4135714817, 3588966466, 457288218, 2540477642, 3898673586,
                                    3750381887, 763821789, 87067
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3073835009, 1333393069, 3746045378, 19326416
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3214557185, 288086972]),
                                        left: Box::new(ProductTree::Leaf(1112338433)),
                                        right: Box::new(ProductTree::Leaf(1112363009))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1201455105, 288129406]),
                                        left: Box::new(ProductTree::Leaf(1112371201)),
                                        right: Box::new(ProductTree::Leaf(1112494081))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1061879809, 3534625656, 3265046213, 19349198
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2040700929, 288229139]),
                                        left: Box::new(ProductTree::Leaf(1112567809)),
                                        right: Box::new(ProductTree::Leaf(1112682497))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3383255041, 288326766]),
                                        left: Box::new(ProductTree::Leaf(1112739841)),
                                        right: Box::new(ProductTree::Leaf(1112887297))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    3336617985, 1472185625, 429406355, 4238434958, 1366729879,
                                    2130923257, 1303727988, 87554
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        996229121, 3895166720, 1095221824, 19374567
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3148677121, 288405306]),
                                        left: Box::new(ProductTree::Leaf(1112952833)),
                                        right: Box::new(ProductTree::Leaf(1112977409))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1068777473, 288528438]),
                                        left: Box::new(ProductTree::Leaf(1113059329)),
                                        right: Box::new(ProductTree::Leaf(1113346049))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        729776129, 679007594, 3364365793, 19409097
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1706876929, 288674965]),
                                        left: Box::new(ProductTree::Leaf(1113403393)),
                                        right: Box::new(ProductTree::Leaf(1113567233))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[4190281729, 288772668]),
                                        left: Box::new(ProductTree::Leaf(1113600001)),
                                        right: Box::new(ProductTree::Leaf(1113747457))
                                    }))
                                }))
                            }))
                        }))
                    })),
                    right: Box::new(ProductTree::Branch(ProductBranch {
                        product: BigUint::from_slice(&[
                            3748446209, 2990918419, 690669137, 4081075455, 1081568834, 933064507,
                            2123731392, 3870738249, 3558696270, 893550403, 3367696824, 1886426235,
                            2514237099, 3695762257, 3253224623, 3917647708, 2484809480, 60790812,
                            841830080, 720782960, 2840955097, 3010593806, 1083663715, 2407562350,
                            672967545, 3363119884, 200783054, 3717269961, 3542047257, 1627244561,
                            3
                        ]),
                        left: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                3942268929, 1068958742, 239300688, 3392926053, 1104372839,
                                3470206431, 249754342, 2190778827, 2595258779, 832880076,
                                3184641538, 4165444005, 2294278113, 3041447333, 3513129662, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    4215619585, 3781710718, 2593826823, 1741148356, 2343840246,
                                    3128730632, 2636066280, 88073
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3314900993, 2186579606, 403551567, 19434097
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3888726017, 288885261]),
                                        left: Box::new(ProductTree::Leaf(1113886721)),
                                        right: Box::new(ProductTree::Leaf(1113894913))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3050053633, 288934128]),
                                        left: Box::new(ProductTree::Leaf(1113935873)),
                                        right: Box::new(ProductTree::Leaf(1114034177))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3786399745, 1583054216, 3987325858, 19464412
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[399687681, 289046752]),
                                        left: Box::new(ProductTree::Leaf(1114165249)),
                                        right: Box::new(ProductTree::Leaf(1114238977))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2984058881, 289223166]),
                                        left: Box::new(ProductTree::Leaf(1114451969)),
                                        right: Box::new(ProductTree::Leaf(1114632193))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    2947874817, 3735384545, 424749002, 2544201653, 772295215,
                                    3938938913, 917023197, 88654
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        3989676033, 3572157506, 3096234457, 19498490
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[4293115905, 289335847]),
                                        left: Box::new(ProductTree::Leaf(1114746881)),
                                        right: Box::new(ProductTree::Leaf(1114771457))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3857309697, 289440041]),
                                        left: Box::new(ProductTree::Leaf(1114894337)),
                                        right: Box::new(ProductTree::Leaf(1115025409))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1172791297, 1740295110, 2514793488, 19528021
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3287187457, 289518732]),
                                        left: Box::new(ProductTree::Leaf(1115074561)),
                                        right: Box::new(ProductTree::Leaf(1115148289))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[3522748417, 289695292]),
                                        left: Box::new(ProductTree::Leaf(1115410433)),
                                        right: Box::new(ProductTree::Leaf(1115492353))
                                    }))
                                }))
                            }))
                        })),
                        right: Box::new(ProductTree::Branch(ProductBranch {
                            product: BigUint::from_slice(&[
                                1416790017, 105566479, 3834809418, 3685214468, 2301695708,
                                2451682437, 443202805, 2842765284, 4273768988, 2248467056,
                                57125895, 2433687461, 3126043589, 3423539492, 3687668514, 1
                            ]),
                            left: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    2114199553, 1475856848, 1980947610, 4155264153, 343470500,
                                    2768905601, 1494286905, 89067
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        637345793, 1026724691, 2675355013, 19552992
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[201162753, 289774018]),
                                        left: Box::new(ProductTree::Leaf(1115590657)),
                                        right: Box::new(ProductTree::Leaf(1115615233))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[704618497, 289810192]),
                                        left: Box::new(ProductTree::Leaf(1115631617)),
                                        right: Box::new(ProductTree::Leaf(1115713537))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        2013724673, 2544985489, 2066258859, 19564337
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1442963457, 289848497]),
                                        left: Box::new(ProductTree::Leaf(1115729921)),
                                        right: Box::new(ProductTree::Leaf(1115762689))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2248482817, 289903830]),
                                        left: Box::new(ProductTree::Leaf(1115803649)),
                                        right: Box::new(ProductTree::Leaf(1115901953))
                                    }))
                                }))
                            })),
                            right: Box::new(ProductTree::Branch(ProductBranch {
                                product: BigUint::from_slice(&[
                                    2322489345, 4263876855, 2781568538, 2933266323, 395075444,
                                    369396914, 3077558088, 89624
                                ]),
                                left: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1579687937, 2831389217, 2800447846, 19602428
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[1745862657, 290084765]),
                                        left: Box::new(ProductTree::Leaf(1116131329)),
                                        right: Box::new(ProductTree::Leaf(1116270593))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[236478465, 290231684]),
                                        left: Box::new(ProductTree::Leaf(1116418049)),
                                        right: Box::new(ProductTree::Leaf(1116549121))
                                    }))
                                })),
                                right: Box::new(ProductTree::Branch(ProductBranch {
                                    product: BigUint::from_slice(&[
                                        1682325505, 1669325779, 518370734, 19637119
                                    ]),
                                    left: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[740261889, 290353079]),
                                        left: Box::new(ProductTree::Leaf(1116639233)),
                                        right: Box::new(ProductTree::Leaf(1116794881))
                                    })),
                                    right: Box::new(ProductTree::Branch(ProductBranch {
                                        product: BigUint::from_slice(&[2552676353, 290476631]),
                                        left: Box::new(ProductTree::Leaf(1116917761)),
                                        right: Box::new(ProductTree::Leaf(1116991489))
                                    }))
                                }))
                            }))
                        }))
                    }))
                }))
            }))
        }))
    });
}

#[cfg(test)]
mod test {

    use std::ops::Mul;

    use itertools::Itertools;
    use num::BigUint;
    use num::FromPrimitive;
    use num::One;
    use num::ToPrimitive;
    use proptest::arbitrary::Arbitrary;
    use proptest::collection::vec;
    use proptest::prop_assert_eq;
    use proptest::strategy::BoxedStrategy;
    use proptest::strategy::Just;
    use proptest::strategy::Strategy;
    use test_strategy::proptest as strategy_proptest;

    use crate::product_tree::ProductTree;
    use crate::residue_number_system::MODULI;

    use super::MASTER_TREE;

    fn arbitrary_biguint(bitlen: usize) -> BoxedStrategy<BigUint> {
        let limbs = vec(u32::arbitrary(), bitlen.div_ceil(32));
        limbs
            .prop_map(move |bb| BigUint::from_slice(&bb) >> (bb.len() * 32 - bitlen))
            .boxed()
    }

    #[strategy_proptest]
    fn cumulative_product(
        #[strategy(1usize..5)] _logn: usize,
        #[strategy(Just(1<<#_logn))] _n: usize,
        #[strategy(vec(u32::arbitrary(), #_n))] ints: Vec<u32>,
    ) {
        let cumulative_product = ints
            .iter()
            .cloned()
            .map(|i| BigUint::from_u32(i).unwrap())
            .fold(BigUint::one(), BigUint::mul);
        let tree = ProductTree::from_leafs(&ints);
        prop_assert_eq!(cumulative_product, tree.value());
    }

    #[strategy_proptest]
    fn reduce_individually(
        #[strategy(1usize..5)] _logn: usize,
        #[strategy(Just(1<<#_logn))] _n: usize,
        #[strategy(vec(u32::arbitrary(), #_n))] moduli: Vec<u32>,
        #[strategy(arbitrary_biguint(#_n*30))] integer: BigUint,
    ) {
        let tree = ProductTree::from_leafs(&moduli);
        let leaf_remainders = tree.reduce(&integer);
        let individual_remainders = moduli
            .into_iter()
            .map(|m| BigUint::from_u32(m).unwrap())
            .map(|p| integer.clone() % p)
            .map(|b| b.to_u32().unwrap())
            .collect_vec();
        prop_assert_eq!(individual_remainders, leaf_remainders);
    }

    #[test]
    fn static_product_tree_matches_with_moduli() {
        let fresh_tree = ProductTree::from_leafs(&MODULI);
        assert_eq!(*MASTER_TREE, fresh_tree);
    }
}
