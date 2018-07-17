#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <mutex>
#include <functional>
#include "seal.h"
#include <random>
#include "util/exring.h"
#include "util/expolycrt.h"
#include <thread>

#ifdef __GNUG__
#include <x86intrin.h>
#endif

using namespace std;
using namespace seal;
using namespace util;

void print_example_banner(string title);
void example_timing();
void example_basics();
void example_weighted_average();
void example_parameter_selection();
void example_batching();
void example_relinearization();
void example_relinearization_part1();
void example_relinearization_part2();
void example_timing();
void example_mt_timing(int th_count);
void example_batching_timing();
void rns_variant_test();
void fv_test_with_rns();
void example_timing_rns();
void example_applygalois();

void test_for_added_functions();
void test_for_recryptor();
void example_floor();
void example_sigmoid();

void test_centereddigits(); 
void test_digitextraction(); 
void test_mulbylearningrate(); 

void test_combineexpand();

int main()
{
	//test_for_added_functions();
	test_for_recryptor();
	//test_centereddigits();
	//example_floor();
	//test_digitextraction();
	//example_sigmoid();
	//example_applygalois();
	//test_combineexpand();
	//test_mulbylearningrate();
	//rns_variant_test();
    // Example: RNS vs FV timings
 //   example_timing_rns();

    //example_timing();

    //rns_variant_test();
    /*
    // Example: Basics
    example_basics();

    // Example: Weighted Average
    example_weighted_average();

    // Example: Automatic Parameter Selection
    example_parameter_selection();

    // Example: Batching using CRT
    example_batching();

    // Example: Relinearization
    example_relinearization();

    // Example: Timing of basic operations
    example_timing();

    // Example: Multi-threaded timing of basic operations
    example_mt_timing(2);
    */
    // Example: Timing of different batching methods
    //example_batching_timing();

    // Wait for ENTER before closing screen.
    //cout << "Press ENTER to exit" << endl;
    //char ignore;
    //cin.get(ignore);

    return 0;
}

void example_basics()
{
    print_example_banner("Example: Basics");

    /*
    In this example we demonstrate using some of the basic arithmetic operations on integers.

    SEAL uses the Fan-Vercauteren (FV) homomorphic encryption scheme. We refer to
    https://eprint.iacr.org/2012/144 for full details on how the FV scheme works.
    */

    // Create encryption parameters.
    EncryptionParameters parms;

    /*
    We first choose the polynomial modulus. This must be a power-of-2 cyclotomic polynomial,
    i.e. a polynomial of the form "1x^(power-of-2) + 1". We recommend using polynomials of
    degree at least 2048.
    */
    parms.set_poly_modulus("1x^2048 + 1");

    /*
    Next we choose the coefficient modulus. SEAL comes with default values for the coefficient
    modulus for some of the most reasonable choices of poly_modulus. They are as follows:

    /----------------------------------------------------------------------\
    | poly_modulus | default coeff_modulus                      | security |
    | -------------|--------------------------------------------|----------|
    | 1x^2048 + 1  | 2^60 - 2^14 + 1 (60 bits)                  | 119 bit  |
    | 1x^4096 + 1  | 2^116 - 2^18 + 1 (116 bits)                | 122 bit  |
    | 1x^8192 + 1  | 2^226 - 2^26 + 1 (226 bits)                | 124 bit  |
    | 1x^16384 + 1 | 2^435 - 2^33 + 1 (435 bits)                | 130 bit  |
    | 1x^32768 + 1 | 2^889 - 2^54 - 2^53 - 2^52 + 1 (889 bits)  | 127 bit  |
    \----------------------------------------------------------------------/

    These can be conveniently accessed using ChooserEvaluator::default_parameter_options(), 
    which returns the above list of options as an std::map, keyed by the degree of the polynomial 
    modulus. The security levels are estimated based on https://eprint.iacr.org/2015/046 and
    https://eprint.iacr.org/2017/047. We strongly recommend that the user consult an expert in 
    the security of RLWE-based cryptography to estimate the security of a particular choice
    of parameters.
    
    The user can also easily choose their custom coefficient modulus. For best performance, it should 
    be a prime of the form 2^A - B, where B is congruent to 1 modulo 2*degree(poly_modulus), and as small 
    as possible. Roughly speaking, When the rest of the parameters are held fixed, increasing coeff_modulus
    decreases the security level. Thus we would not recommend using a value for coeff_modulus much larger 
    than those listed above (the defaults). In general, we highly recommend the user to consult with an expert 
    in the security of RLWE-based cryptography when selecting their parameters to ensure an appropriate level 
    of security.

    The size of coeff_modulus affects the total noise budget that a freshly encrypted ciphertext has. More 
    precisely, every ciphertext starts with a certain amount of noise budget, which is consumed in homomorphic
    operations - in particular in multiplication. Once the noise budget reaches 0, the ciphertext becomes 
    impossible to decrypt. The total noise budget in a freshly encrypted ciphertext is very roughly given by 
    log2(coeff_modulus/plain_modulus), so increasing coeff_modulus will allow the user to perform more
    homomorphic operations on the ciphertexts without corrupting them. However, we must again warn that
    increasing coeff_modulus has a strong negative effect on the security level.
    */
    parms.set_coeff_modulus(ChooserEvaluator::default_parameter_options().at(2048));

    /*
    Now we set the plaintext modulus. This can be any positive integer, even though here we take it to be a 
    power of two. A larger plaintext modulus causes the noise to grow faster in homomorphic multiplication, 
    and also lowers the maximum amount of noise in ciphertexts that the system can tolerate (see above).
    On the other hand, a larger plaintext modulus typically allows for better homomorphic integer arithmetic,
    although this depends strongly on which encoder is used to encode integers into plaintext polynomials.
    */
    parms.set_plain_modulus(1 << 8);

    /*
    Once all parameters are set, we need to create a SEALContext object, which evaluates the properties 
    of the parameters, their validity for homomorphic encryption, and performs some important pre-computation.
    */
    SEALContext context(parms);

    /*
    Plaintext elements in the FV scheme are polynomials (represented by the Plaintext class) with coefficients 
    integers modulo plain_modulus. To encrypt for example integers instead, one must use an "encoding scheme", 
    i.e. a specific way of representing integers as such polynomials. SEAL comes with a few basic encoders:

    IntegerEncoder:
    Given an integer base b, encodes integers as plaintext polynomials in the following way. First, a base-b
    expansion of the integer is computed. This expansion uses a "balanced" set of representatives of integers
    modulo b as the coefficients. Namely, when b is off the coefficients are integers between -(b-1)/2 and
    (b-1)/2. When b is even, the integers are between -b/2 and (b-1)/2, except when b is two and the usual
    binary expansion is used (coefficients 0 and 1). Decoding amounts to evaluating the polynomial at x=b.
    For example, if b=2, the integer 26 = 2^4 + 2^3 + 2^1 is encoded as the polynomial 1x^4 + 1x^3 + 1x^1.
    When b=3, 26 = 3^3 - 3^0 is encoded as the polynomial 1x^3 - 1. In reality, coefficients of polynomials
    are always unsigned integers, and in this case are stored as their smallest non-negative representatives
    modulo plain_modulus. To create an integer encoder with a base b, use IntegerEncoder(plain_modulus, b). 
    If no b is given to the constructor, the default value of b=2 is used.

    FractionalEncoder:
    Encodes fixed-precision rational numbers as follows. First expand the number in a given base b, possibly 
    truncating an infinite fractional part to finite precision, e.g. 26.75 = 2^4 + 2^3 + 2^1 + 2^(-1) + 2^(-2)
    when b=2. For the sake of the example, suppose poly_modulus is 1x^1024 + 1. Next represent the integer part 
    of the number in the same way as in IntegerEncoder (with b=2 here). Finally, represent the fractional part 
    in the leading coefficients of the polynomial, but when doing so invert the signs of the coefficients. So 
    in this example we would represent 26.75 as the polynomial -1x^1023 - 1x^1022 + 1x^4 + 1x^3 + 1x^1. The 
    negative coefficients of the polynomial will again be represented as their negatives modulo plain_modulus.

    PolyCRTBuilder:
    If poly_modulus is 1x^N + 1, PolyCRTBuilder allows "batching" of N plaintext integers modulo plain_modulus 
    into one plaintext polynomial, where homomorphic operations can be carried out very efficiently in a SIMD 
    manner by operating on such a "composed" plaintext or ciphertext polynomials. For full details on this very
    powerful technique we recommend https://eprint.iacr.org/2012/565.pdf and https://eprint.iacr.org/2011/133.

    A crucial fact to understand is that when homomorphic operations are performed on ciphertexts, they will
    carry over to the underlying plaintexts, and as a result of additions and multiplications the coefficients
    in the plaintext polynomials will increase from what they originally were in freshly encoded polynomials.
    This becomes a problem when the coefficients reach the size of plain_modulus, in which case they will get
    automatically reduced modulo plain_modulus, and might render the underlying plaintext polynomial impossible
    to be correctly decoded back into an integer or rational number. Therefore, it is typically crucial to
    have a good sense of how large the coefficients will grow in the underlying plaintext polynomials when
    homomorphic computations are carried out on the ciphertexts, and make sure that plain_modulus is chosen to
    be at least as large as this number.

    Here we choose to create an IntegerEncoder with base b=2.
    */
    IntegerEncoder encoder(context.plain_modulus());

    // Encode two integers as polynomials.
    const int value1 = 5;
    const int value2 = -7;
    Plaintext encoded1 = encoder.encode(value1);
    Plaintext encoded2 = encoder.encode(value2);
    cout << "Encoded " << value1 << " as polynomial " << encoded1.to_string() << endl;
    cout << "Encoded " << value2 << " as polynomial " << encoded2.to_string() << endl;
    
    // Generate keys.
    cout << "Generating keys ..." << endl;
    KeyGenerator generator(context);
    generator.generate();
    cout << "... key generation complete" << endl;
    PublicKey public_key = generator.public_key();
    SecretKey secret_key = generator.secret_key();

    // Encrypt values.
    cout << "Encrypting values..." << endl;
    Encryptor encryptor(context, public_key);
    Ciphertext encrypted1 = encryptor.encrypt(encoded1);
    Ciphertext encrypted2 = encryptor.encrypt(encoded2);

    // Perform arithmetic on encrypted values.
    cout << "Performing arithmetic on ecrypted numbers ..." << endl;
    Evaluator evaluator(context);
    cout << "Performing homomorphic negation ..." << endl;
    Ciphertext encryptednegated1 = evaluator.negate(encrypted1);
    cout << "Performing homomorphic addition ..." << endl;
    Ciphertext encryptedsum = evaluator.add(encrypted1, encrypted2);
    cout << "Performing homomorphic subtraction ..." << endl;
    Ciphertext encrypteddiff = evaluator.sub(encrypted1, encrypted2);
    cout << "Performing homomorphic multiplication ..." << endl;
    Ciphertext encryptedproduct = evaluator.multiply(encrypted1, encrypted2);

    // Decrypt results.
    cout << "Decrypting results ..." << endl;
    Decryptor decryptor(context, secret_key);
    Plaintext decrypted1 = decryptor.decrypt(encrypted1);
    Plaintext decrypted2 = decryptor.decrypt(encrypted2);
    Plaintext decryptednegated1 = decryptor.decrypt(encryptednegated1);
    Plaintext decryptedsum = decryptor.decrypt(encryptedsum);
    Plaintext decrypteddiff = decryptor.decrypt(encrypteddiff);
    Plaintext decryptedproduct = decryptor.decrypt(encryptedproduct);

    // Decode results.
    int decoded1 = encoder.decode_int32(decrypted1);
    int decoded2 = encoder.decode_int32(decrypted2);
    int decodednegated1 = encoder.decode_int32(decryptednegated1);
    int decodedsum = encoder.decode_int32(decryptedsum);
    int decodeddiff = encoder.decode_int32(decrypteddiff);
    int decodedproduct = encoder.decode_int32(decryptedproduct);

    // Display results.
    cout << "Original = " << value1 << "; after encryption/decryption = " << decoded1 << endl;
    cout << "Original = " << value2 << "; after encryption/decryption = " << decoded2 << endl;
    cout << "Encrypted negate of " << value1 << " = " << decodednegated1 << endl;
    cout << "Encrypted addition of " << value1 << " and " << value2 << " = " << decodedsum << endl;
    cout << "Encrypted subtraction of " << value1 << " and " << value2 << " = " << decodeddiff << endl;
    cout << "Encrypted multiplication of " << value1 << " and " << value2 << " = " << decodedproduct << endl;

    // How much noise budget did we use in these operations?
    cout << "Noise budget in encryption of " << value1 << ": " 
        << decryptor.invariant_noise_budget(encrypted1) << " bits" << endl;
    cout << "Noise budget in encryption of " << value2 << ": "
        << decryptor.invariant_noise_budget(encrypted2) << " bits" << endl;
    cout << "Noise budget in sum: "
        << decryptor.invariant_noise_budget(encryptedsum) << " bits" << endl;
    cout << "Noise budget in product: "
        << decryptor.invariant_noise_budget(encryptedproduct) << " bits" << endl;
}

void example_weighted_average()
{
    print_example_banner("Example: Weighted Average");

    // In this example we demonstrate computing a weighted average of 10 rational numbers.
    
    // The 10 rational numbers we use are:
    const vector<double> rational_numbers { 3.1, 4.159, 2.65, 3.5897, 9.3, 2.3, 8.46, 2.64, 3.383, 2.795 };

    // The 10 weights are:
    const vector<double> coefficients { 0.1, 0.05, 0.05, 0.2, 0.05, 0.3, 0.1, 0.025, 0.075, 0.05 };

    // Create encryption parameters
    EncryptionParameters parms;
    parms.set_poly_modulus("1x^2048 + 1");
    parms.set_coeff_modulus(ChooserEvaluator::default_parameter_options().at(2048));
    parms.set_plain_modulus(1 << 8);
    SEALContext context(parms);
    // Generate keys.
    cout << "Generating keys ..." << endl;
    KeyGenerator generator(context);
    generator.generate();
    cout << "... key generation complete" << endl;
    PublicKey public_key = generator.public_key();
    SecretKey secret_key = generator.secret_key();

    /*
    We will need a fractional encoder for dealing with the rational numbers. Here we reserve 
    64 coefficients of the polynomial for the integral part (low-degree terms) and expand the 
    fractional part to 32 terms of precision (base 3) (high-degree terms).
    */
    FractionalEncoder encoder(context.plain_modulus(), context.poly_modulus(), 64, 32, 3);

    // Create the rest of the tools
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    // First we encrypt the rational numbers
    cout << "Encrypting ... ";
    vector<Ciphertext> encrypted_rationals;
    for (int i = 0; i < 10; i++)
    {
        Plaintext encoded_number = encoder.encode(rational_numbers[i]);
        encrypted_rationals.emplace_back(encryptor.encrypt(encoded_number));
        cout << to_string(rational_numbers[i]).substr(0,6) << ((i < 9) ? ", " : ".\n");
    }

    // Next we encode the coefficients. There is no reason to encrypt these since they are not private data.
    cout << "Encoding ... ";
    vector<Plaintext> encoded_coefficients;
    for (int i = 0; i < 10; i++)
    {
        encoded_coefficients.emplace_back(encoder.encode(coefficients[i]));
        cout << to_string(coefficients[i]).substr(0,6) << ((i < 9) ? ", " : ".\n");
    }
    
    // We also need to encode 0.1. We will multiply the result by this to perform division by 10.
    Plaintext div_by_ten = encoder.encode(0.1);

    // Now compute all the products of the encrypted rational numbers with the plaintext coefficients
    cout << "Computing products ... ";
    vector<Ciphertext> encrypted_products;
    for (int i = 0; i < 10; i++)
    {
        /*
        We use Evaluator::multiply_plain(...) instead of Evaluator::multiply(...) (which would 
        require also the coefficient to be encrypted). This has much better noise growth
        behavior than multiplying two encrypted numbers does.
        */
        Ciphertext enc_plain_product = evaluator.multiply_plain(encrypted_rationals[i], encoded_coefficients[i]);
        encrypted_products.emplace_back(enc_plain_product);
    }
    cout << "done." << endl;

    // Now we add together these products. The most convenient way to do that is
    // to use the function Evaluator::add_many(...).
    cout << "Add up all 10 ciphertexts ... ";
    Ciphertext encrypted_dot_product = evaluator.add_many(encrypted_products);
    cout << " done." << endl;

    // Finally we divide by 10 to obtain the result.
    cout << "Divide by 10 ... ";
    Ciphertext encrypted_result = evaluator.multiply_plain(encrypted_dot_product, div_by_ten);
    cout << "done." << endl;

    // Decrypt
    cout << "Decrypting ... ";
    Plaintext plain_result(decryptor.decrypt(encrypted_result));
    cout << "done." << endl;

    // Print the result
    double result = encoder.decode(plain_result);
    cout << "Weighted average: " << result << endl;

    // How much noise budget are we left with?
    cout << "Noise budget in result: " << decryptor.invariant_noise_budget(encrypted_result) << " bits" << endl;
}

void example_parameter_selection()
{
    print_example_banner("Example: Automatic Parameter Selection");

    /*
    Here we demonstrate the automatic parameter selection tool. Suppose we want to find parameters
    that are optimized in a way that allows us to evaluate the polynomial 42x^3-27x+1. We need to know
    the size of the input data, so let's assume that x is an integer with base-3 representation of length
    at most 10.
    */
    cout << "Finding optimized parameters for computing 42x^3-27x+1 ... ";

    ChooserEncoder chooser_encoder(3);
    ChooserEvaluator chooser_evaluator;
    
    /*
    First create a ChooserPoly representing the input data. You can think of this modeling a freshly
    encrypted cipheretext of a plaintext polynomial with length at most 10 coefficients, where the
    coefficients have absolute value at most 1.
    */
    ChooserPoly cinput(10, 1);

    // Compute the first term
    ChooserPoly ccubed_input = chooser_evaluator.exponentiate(cinput, 3);
    ChooserPoly cterm1 = chooser_evaluator.multiply_plain(ccubed_input, chooser_encoder.encode(42));

    // Compute the second term
    ChooserPoly cterm2 = chooser_evaluator.multiply_plain(cinput, chooser_encoder.encode(27));

    // Subtract the first two terms
    ChooserPoly csum12 = chooser_evaluator.sub(cterm1, cterm2);

    // Add the constant term 1
    ChooserPoly cresult = chooser_evaluator.add_plain(csum12, chooser_encoder.encode(1));

    // To find an optimized set of parameters, we use ChooserEvaluator::select_parameters(...).
    EncryptionParameters optimal_parms;
    chooser_evaluator.select_parameters({ cresult }, 0, optimal_parms);
    
    // Create an SEALContext object for the returned parameters
    SEALContext optimal_context(optimal_parms);

    cout << "done." << endl;

    // Let's print these to see what was recommended
    cout << "Selected parameters:" << endl;
    cout << "| poly_modulus: " << optimal_context.poly_modulus().to_string() << endl;
    cout << "| coeff_modulus: " << optimal_context.coeff_modulus().to_string() << endl;
    cout << "| plain_modulus: " << optimal_context.plain_modulus().to_dec_string() << endl;
    cout << "| decomposition_bit_count: " << optimal_context.decomposition_bit_count() << endl;
    cout << "| noise_standard_deviation: " << optimal_context.noise_standard_deviation() << endl;
    cout << "| noise_max_deviation: " << optimal_context.noise_max_deviation() << endl;

    // Let's try to actually perform the homomorphic computation using the recommended parameters.
    // Generate keys.
    cout << "Generating keys ..." << endl;
    KeyGenerator generator(optimal_context);

    /*
    Need to generate one evaluation key because below we will use Evaluator::exponentiate(...),
    which relinearizes after every multiplication it performs (see example_relinearization()
    for more details.
    */
    generator.generate(1);
    cout << "... key generation complete" << endl;
    PublicKey public_key = generator.public_key();
    SecretKey secret_key = generator.secret_key();
    EvaluationKeys evk = generator.evaluation_keys();

    // Create the encoding/encryption tools
    IntegerEncoder encoder(optimal_context.plain_modulus(), 3);
    Encryptor encryptor(optimal_context, public_key);
    Evaluator evaluator(optimal_context, evk);
    Decryptor decryptor(optimal_context, secret_key);

    // Now perform the computations on real encrypted data.
    int input_value = 12345;
    Plaintext plain_input = encoder.encode(input_value);
    cout << "Encoded " << input_value << " as polynomial " << plain_input.to_string() << endl;

    cout << "Encrypting ... ";
    Ciphertext input = encryptor.encrypt(plain_input);
    cout << "done." << endl;

    // Compute the first term
    cout << "Computing first term ... ";
    Ciphertext cubed_input = evaluator.exponentiate(input, 3);
    Ciphertext term1 = evaluator.multiply_plain(cubed_input, encoder.encode(42));
    cout << "done." << endl;

    // Compute the second term
    cout << "Computing second term ... ";
    Ciphertext term2 = evaluator.multiply_plain(input, encoder.encode(27));
    cout << "done." << endl;

    // Subtract the first two terms
    cout << "Subtracting first two terms ... ";
    Ciphertext sum12 = evaluator.sub(term1, term2);
    cout << "done." << endl;

    // Add the constant term 1
    cout << "Adding one ... ";
    Ciphertext result = evaluator.add_plain(sum12, encoder.encode(1));
    cout << "done." << endl;

    // Decrypt and decode
    cout << "Decrypting ... ";
    Plaintext plain_result = decryptor.decrypt(result);
    cout << "done." << endl;
    
    // Finally print the result
    cout << "Polynomial 42x^3-27x+1 evaluated at x=12345: " << encoder.decode_int64(plain_result) << endl;

    // How much noise budget are we left with?
    cout << "Noise budget in result: " << decryptor.invariant_noise_budget(result) << " bits" << endl;
}

void example_batching()
{
    print_example_banner("Example: Batching using CRT");

    // Create encryption parameters
    EncryptionParameters parms;

    /*
    For PolyCRTBuilder it is necessary to have plain_modulus be a prime number congruent to 1 modulo 
    2*degree(poly_modulus). We can use for example the following parameters:
    */
    parms.set_poly_modulus("1x^4096 + 1");
    parms.set_coeff_modulus(ChooserEvaluator::default_parameter_options().at(4096));
    parms.set_plain_modulus(40961);
    SEALContext context(parms);

    // Create the PolyCRTBuilder
    PolyCRTBuilder crtbuilder(context);
    int slot_count = crtbuilder.get_slot_count();

    cout << "Encryption parameters allow " << slot_count << " slots." << endl;

    // Create a vector of values that are to be stored in the slots. We initialize all values to 0 at this point.
    vector<BigUInt> values(slot_count, BigUInt(context.plain_modulus().bit_count(), static_cast<uint64_t>(0)));

    // Set the first few entries of the values vector to be non-zero
    values[0] = 2;
    values[1] = 3;
    values[2] = 5;
    values[3] = 7;
    values[4] = 11;
    values[5] = 13;

    // Now compose these into one polynomial using PolyCRTBuilder
    cout << "Plaintext slot contents (slot, value): ";
    for (size_t i = 0; i < 6; i++)
    {
        cout << "(" << i << ", " << values[i].to_dec_string() << ")" << ((i != 5) ? ", " : "\n");
    }
    Plaintext plain_composed_poly = crtbuilder.compose(values);

    // Let's do some homomorphic operations now. First we need all the encryption tools.
    // Generate keys.
    cout << "Generating keys ..." << endl;
    KeyGenerator generator(context);
    generator.generate();
    cout << "... key generation complete" << endl;
    PublicKey public_key = generator.public_key();
    SecretKey secret_key = generator.secret_key();

    // Create the encryption tools
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);

    // Encrypt plain_composed_poly
    cout << "Encrypting ... ";
    Ciphertext encrypted_composed_poly = encryptor.encrypt(plain_composed_poly);
    cout << "done." << endl;

    // Let's square the encrypted_composed_poly
    cout << "Squaring the encrypted polynomial ... ";
    Ciphertext encrypted_square = evaluator.square(encrypted_composed_poly);
    cout << "done." << endl;

    cout << "Decrypting the squared polynomial ... ";
    Plaintext plain_square = decryptor.decrypt(encrypted_square);
    cout << "done." << endl;
    
    // Print the squared slots
    crtbuilder.decompose(plain_square, values);
    cout << "Squared slot contents (slot, value): ";
    for (size_t i = 0; i < 6; i++)
    {
        cout << "(" << i << ", " << values[i].to_dec_string() << ")" << ((i != 5) ? ", " : "\n");
    }

    // Now let's try to multiply the squares with the plaintext coefficients (3, 1, 4, 1, 5, 9, 0, 0, ..., 0).
    // First create the coefficient vector
    vector<BigUInt> plain_coeff_vector(slot_count, BigUInt(context.plain_modulus().bit_count(), static_cast<uint64_t>(0)));
    plain_coeff_vector[0] = 3;
    plain_coeff_vector[1] = 1;
    plain_coeff_vector[2] = 4;
    plain_coeff_vector[3] = 1;
    plain_coeff_vector[4] = 5;
    plain_coeff_vector[5] = 9;

    // Use PolyCRTBuilder to compose plain_coeff_vector into a polynomial
    Plaintext plain_coeff_poly = crtbuilder.compose(plain_coeff_vector);

    // Print the coefficient vector
    cout << "Coefficient slot contents (slot, value): ";
    for (size_t i = 0; i < 6; i++)
    {
        cout << "(" << i << ", " << plain_coeff_vector[i].to_dec_string() << ")" << ((i != 5) ? ", " : "\n");
    }

    // Now use multiply_plain to multiply each encrypted slot with the corresponding coefficient
    cout << "Multiplying squared slots with the coefficients ... ";
    Ciphertext encrypted_scaled_square = evaluator.multiply_plain(encrypted_square, plain_coeff_poly);
    cout << " done." << endl;
    
    // Decrypt it
    cout << "Decrypting the scaled squared polynomial ... ";
    Plaintext plain_scaled_square = decryptor.decrypt(encrypted_scaled_square);
    cout << "done." << endl;

    // Print the scaled squared slots
    cout << "Scaled squared slot contents (slot, value): ";
    crtbuilder.decompose(plain_scaled_square, values);
    for (size_t i = 0; i < 6; i++)
    {
        cout << "(" << i << ", " << values[i].to_dec_string() << ")" << ((i != 5) ? ", " : "\n");
    }

    // How much noise budget are we left with?
    cout << "Noise budget in result: " << decryptor.invariant_noise_budget(encrypted_scaled_square) << " bits" << endl;
}

void example_relinearization()
{
    print_example_banner("Example: Relinearization");

    /*
    A valid ciphertext consists of at least two polynomials. To read the current size of a ciphertext 
    the user can use Ciphertext::size(). A fresh ciphertext always has size 2, and performing 
    homomorphic multiplication results in the output ciphertext growing in size. More precisely, 
    if the input ciphertexts have size M and N, then the output ciphertext after homomorphic 
    multiplication will have size M+N-1.

    The multiplication operation on input ciphertexts of size M and N will require O(M*N) polynomial 
    multiplications to be performed. Therefore, the multiplication of large ciphertexts could be 
    very computationally costly and in some situations the user might prefer to reduce the size of 
    the ciphertexts by performing a so-called relinearization operation.

    The function Evaluator::relinearize(...) can reduce the size of an input ciphertext of size M 
    to any size in 2, 3, ..., M-1. As was explained above, relinearizing one or both of two ciphertexts 
    before performing multiplication on them may significantly reduce the computational cost of the 
    multiplication. However, note that the relinearization process itself also requires several polynomial 
    multiplications to be performed. Using the Number Theoretic Transform (NTT) for relinearization, 
    reducing a ciphertext of size K to size L requires (K - L)*([floor(log_2(coeff_modulus)/dbc) + 3]) 
    NTT transforms, where dbc denotes the encryption parameter "decomposition_bit_count".

    Relinearization also affects the inherent noise in two ways. First, a larger ciphertexts produces
    more noise in homomorphic multiplication than a smaller one does. If the ciphertexts are small,
    the effect of the ciphertext size is insignificant, but if they are very large the effect can 
    easily become the biggest contributor to noise. Second, relinearization increases the inherent 
    noise in the ciphertext to be relinearized by an additive factor. This should be contrasted with 
    the multiplicative factor that homomorphic multiplication increases the noise by. The additive
    factor is proportional to 2^dbc, which can be either very small or very large compared to the 
    current level of inherent noise in the ciphertext. This means that if the ciphertextis very fresh 
    (has very little noise in it), relinearization might have a significant adverse effect on the 
    homomorphic computation ability, and it might make sense to instead use larger ciphertexts and 
    relinearize at a later point where the additive noise term vanishes into an already larger noise,
    or alternatively use a smaller dbc, which will result in slightly slower relinearization. 
    
    When using relinearization it is necessary that the decomposition_bit_count variable is set to 
    some positive value in the encryption parameters, and that enough evaluation keys are given to
    the constructor of Evaluator. We will discuss evaluation keys when we construct the key generator.
   
    We will provide two examples of relinearization.
    */

    /*
    Example 1: We demonstrate using Evaluator::relinearize(...) and illustrate how it reduces the 
    ciphertext sizes at the cost of increasing running time and noise in a particular computation.
    */
    example_relinearization_part1();
    cout << endl;

    /*
    Example 2: We demonstrate how relinearization can reduce both running time and noise. 
    */
    example_relinearization_part2();
}

void example_relinearization_part1()
{
    cout << "Example 1: Performing relinearization too early can increase noise in the final result." << endl;

    // Set up encryption parameters
    EncryptionParameters parms;
    parms.set_poly_modulus("1x^4096 + 1");
    parms.set_coeff_modulus(ChooserEvaluator::default_parameter_options().at(4096));
    parms.set_plain_modulus(1 << 8);

    /*
    The choice of decomposition_bit_count (dbc) can affect the performance of relinearization
    noticeably. A reasonable choice for it is between 1/10 and 1/2 of the significant bit count 
    of the coefficient modulus. Sometimes when the dbc needs to be very small (due to noise growth), 
    it might make more sense to move up to a larger poly_modulus and coeff_modulus, and set dbc to 
    be as large as possible.

    A smaller dbc will make relinearization too slow. A higher dbc will increase noise growth 
    while not making relinearization any faster. Here, the coeff_modulus has 116 significant 
    bits, so we choose dbc to be half of this. We can expect to see extreme differences in 
    noise growth between the relinearizing and non-relinearizing cases due to the decomposition 
    bit count being so large.
    */
    parms.set_decomposition_bit_count(58);

    // Done with parameters
    SEALContext context(parms);

    /*
    By default, KeyGenerator::generate() will generate no evaluation keys. This means that we 
    cannot perform any relinearization. However, this is sufficient for performing all other 
    homomorphic evaluation operations as they do not use evaluation keys, and is enough for
    now as we start by demonstrating the computation without relinearization.
    */
    cout << "Generating keys ..." << endl;
    KeyGenerator generator(context);
    generator.generate();
    cout << "... key generation complete" << endl;
    PublicKey public_key = generator.public_key();
    SecretKey secret_key = generator.secret_key();

    /*
    Suppose we want to homomorphically multiply four ciphertexts together. Does it make sense 
    to relinearize at an intermediate stage of the computation?
    */

    // Encrypt plaintexts to generate the four fresh ciphertexts
    Plaintext plain1("5");
    Plaintext plain2("6");
    Plaintext plain3("7");
    Plaintext plain4("8");
    cout << "Encrypting values { 5, 6, 7, 8 } as { encrypted1, encrypted2, encrypted3, encrypted4 }" << endl;
    Encryptor encryptor(context, public_key);
    Ciphertext encrypted1 = encryptor.encrypt(plain1);
    Ciphertext encrypted2 = encryptor.encrypt(plain2);
    Ciphertext encrypted3 = encryptor.encrypt(plain3);
    Ciphertext encrypted4 = encryptor.encrypt(plain4);

    // We need a Decryptor to be able to measure the inherent noise
    Decryptor decryptor(context, secret_key);

    // What are the noise budgets in the four ciphertexts?
    cout << "Noise budgets in the four ciphertexts: "
        << decryptor.invariant_noise_budget(encrypted1) << " bits, "
        << decryptor.invariant_noise_budget(encrypted2) << " bits, "
        << decryptor.invariant_noise_budget(encrypted3) << " bits, "
        << decryptor.invariant_noise_budget(encrypted4) << " bits" << endl;

    // Construct an Evaluator
    Evaluator evaluator(context);

    // Perform first part of computation
    cout << "Computing enc_prod1 as encrypted1*encrypted2 ..." << endl;
    Ciphertext enc_prod1 = evaluator.multiply(encrypted1, encrypted2);
    cout << "Computing enc_prod2 as encrypted3*encrypted4 ..." << endl;
    Ciphertext enc_prod2 = evaluator.multiply(encrypted3, encrypted4);

    // First the result with no relinearization
    cout << endl;
    cout << "Path 1: No relinearization" << endl;

    // Compute product of all four
    cout << "Computing result as enc_prod1*enc_prod2 ..." << endl;
    Ciphertext enc_result = evaluator.multiply(enc_prod1, enc_prod2);

    // Now enc_result has size 5
    cout << "Size of enc_result: " << enc_result.size() << endl;

    // How much noise budget are we left with?
    int noise_budget_norelin = decryptor.invariant_noise_budget(enc_result);
    cout << "Noise budget in enc_result: " << noise_budget_norelin << " bits" << endl;

    /*
    We didn't create any evaluation keys, so we can't relinearize at all with the current
    Evaluator. In general, relinearizing down from size K to any smaller size (but at least 2)
    requires at least K-2 evaluation keys. In this case we wish to relinearize enc_prod1 and
    enc_prod2, which both have size 3. Thus we need only one evaluation key.

    We can create this new evaluation key by calling KeyGenerator::generate_evaluation_keys(...). 
    Alternatively, we could have created it already in the beginning by instead of calling 
    generator.generate(1) instead of generator.generate().

    We will also need a new Evaluator, as the previous one was constructed without enough 
    indeed, any) evaluation keys. It is not possible to add new evaluation keys to a previously 
    created Evaluator.
    */
    generator.generate_evaluation_keys(1);
    EvaluationKeys evaluation_keys = generator.evaluation_keys();
    Evaluator evaluator2(context, evaluation_keys);

    // Now with relinearization
    cout << endl;
    cout << "Path 2: With relinearization" << endl;

    // What if we do intermediate relinearization of enc_prod1 and enc_prod2?
    cout << "Relinearizing enc_prod1 and enc_prod2 to size 2 ..." << endl;
    Ciphertext enc_relin_prod1 = evaluator2.relinearize(enc_prod1);
    Ciphertext enc_relin_prod2 = evaluator2.relinearize(enc_prod2);

    // Now multiply the relinearized products together
    cout << "Computing enc_result as enc_relin_prod1*enc_relin_prod2 ..." << endl;
    enc_result = evaluator2.multiply(enc_relin_prod1, enc_relin_prod2);

    // Now enc_result has size 3
    cout << "Size of enc_result: " << enc_result.size() << endl;

    // How much noise budget are we left with?
    int noise_budget_relin = decryptor.invariant_noise_budget(enc_result);
    cout << "Noise budget in enc_result: " << noise_budget_relin << " bits" << endl;

    /*
    While in this case the noise increased significantly due to relinearization, in other
    computations the situation might be entirely different. Indeed, recall that larger
    ciphertext sizes can have a huge adverse effect on noise growth in multiplication.
    Also recall that homomorphic multiplication is much slower when the ciphertexts are 
    larger.
    */    
}

void example_relinearization_part2()
{
    cout << "Example 2: Effect on running time and noise in computing [(enc1*enc2)*(enc3*enc4)]^2." << endl;

    // Set up encryption parameters
    EncryptionParameters parms;
    parms.set_poly_modulus("1x^4096 + 1");
    parms.set_coeff_modulus(ChooserEvaluator::default_parameter_options().at(4096));
    parms.set_plain_modulus(1 << 6);

    /*
    We use a relatively small decomposition bit count here to avoid significant noise
    growth from the relinearization operation itself. Make this bigger and you will
    see both increased running time and decreased noise.
    */
    parms.set_decomposition_bit_count(16);

    // Done with parameters
    SEALContext context(parms);

    // We generate the encryption keys and one evaluation key.
    cout << "Generating keys ..." << endl;
    KeyGenerator generator(context);
    generator.generate(1);
    cout << "... key generation complete" << endl;
    PublicKey public_key = generator.public_key();
    SecretKey secret_key = generator.secret_key();
    EvaluationKeys evaluation_keys = generator.evaluation_keys();

    // Encrypt plaintexts to generate the four fresh ciphertexts
    Plaintext plain1("4");
    Plaintext plain2("3x^1");
    Plaintext plain3("2x^2");
    Plaintext plain4("1x^3");
    cout << "Encrypting values { 4, 3x^1, 2x^2, x^3 } as { encrypted1, encrypted2, encrypted3, encrypted4 }" << endl;
    Encryptor encryptor(context, public_key);
    Ciphertext encrypted1 = encryptor.encrypt(plain1);
    Ciphertext encrypted2 = encryptor.encrypt(plain2);
    Ciphertext encrypted3 = encryptor.encrypt(plain3);
    Ciphertext encrypted4 = encryptor.encrypt(plain4);

    // We need a Decryptor to be able to measure the inherent noise
    Decryptor decryptor(context, secret_key);

    // Construct an Evaluator
    Evaluator evaluator(context, evaluation_keys);

    cout << "Computing enc_prod12 = encrypted1*encrypted2 ..." << endl;
    Ciphertext enc_prod12 = evaluator.multiply(encrypted1, encrypted2);

    cout << "Computing enc_prod34 = encrypted3*encrypted4 ..." << endl;
    Ciphertext enc_prod34 = evaluator.multiply(encrypted3, encrypted4);

    // First the result with no relinearization
    cout << endl;
    cout << "Path 1: No relinerization" << endl;

    auto time_norelin_start = chrono::high_resolution_clock::now();

    // Compute product of all four
    cout << "Computing enc_prod = enc_prod12*enc_prod34 ..." << endl;
    Ciphertext enc_prod = evaluator.multiply(enc_prod12, enc_prod34);

    cout << "Computing enc_square = [enc_prod]^2 ..." << endl;
    Ciphertext enc_square = evaluator.square(enc_prod);

    auto time_norelin_end = chrono::high_resolution_clock::now();
    cout << "Time (without relinearization): " << chrono::duration_cast<chrono::microseconds>(time_norelin_end - time_norelin_start).count()
        << " microseconds" << endl;

    // Print size and noise budget of result. 
    cout << "Size of enc_square: " << enc_square.size() << endl;
    cout << "Noise budget in enc_square: " << decryptor.invariant_noise_budget(enc_square) << " bits" << endl;

    // Now the same thing but with relinearization
    cout << endl;
    cout << "Path 2: With relinerization" << endl;

    auto time_relin_start = chrono::high_resolution_clock::now();

    cout << "Relinearizing enc_prod12 and enc_prod34 to size 2 ..." << endl;
    Ciphertext enc_relin_prod12 = evaluator.relinearize(enc_prod12);
    Ciphertext enc_relin_prod34 = evaluator.relinearize(enc_prod34);

    // Now multiply the relinearized products together
    cout << "Computing enc_prod = enc_relin_prod12*enc_relin_prod34... " << endl;
    enc_prod = evaluator.multiply(enc_relin_prod12, enc_relin_prod34);

    cout << "Computing enc_square = [enc_prod]^2 ... " << endl;
    enc_square = evaluator.square(enc_prod);

    auto time_relin_end = chrono::high_resolution_clock::now();
    cout << "Time (with relinearization): " << chrono::duration_cast<chrono::microseconds>(time_relin_end - time_relin_start).count()
        << " microseconds" << endl;

    // Print size and noise budget of result. 
    cout << "Size of enc_square: " << enc_square.size() << endl;
    cout << "Noise budget in enc_square: " << decryptor.invariant_noise_budget(enc_square) << " bits" << endl;
}

void example_timing()
{
    print_example_banner("Example: Timing of basic operations");

    auto performance_test = [](SEALContext &context)
    {
        auto poly_modulus = context.poly_modulus();
        auto coeff_modulus = context.coeff_modulus();
        auto plain_modulus = context.plain_modulus();

        KeyGenerator keygen(context);
        keygen.generate(1);
        auto secret_key = keygen.secret_key();
        auto public_key = keygen.public_key();
        auto evk = keygen.evaluation_keys();

        Encryptor encryptor(context, public_key);
        Decryptor decryptor(context, secret_key);
        Evaluator evaluator(context, evk);
        IntegerEncoder encoder(plain_modulus);

        chrono::microseconds time_encode_sum(0);
        chrono::microseconds time_encrypt_sum(0);
        chrono::microseconds time_multiply_sum(0);
        chrono::microseconds time_square_sum(0);
        chrono::microseconds time_relinearize_sum(0);
        chrono::microseconds time_decrypt_sum(0);
        chrono::microseconds time_decode_sum(0);

        cout << "Encryption parameters:" << endl;
        cout << "| poly_modulus: " << poly_modulus.to_string() << endl;
        cout << "| coeff_modulus: " << coeff_modulus.to_string() << endl;
        cout << "| plain_modulus: " << plain_modulus.to_dec_string() << endl;
        cout << "| decomposition_bit_count: " << context.decomposition_bit_count() << endl;
        cout << "| noise_standard_deviation: " << context.noise_standard_deviation() << endl;
        cout << "| noise_max_deviation: " << context.noise_max_deviation() << endl << endl;

        int count = 10;

        cout << "Running tests ";
        for (int i = 0; i < count; i++)
        {
            auto time_start = chrono::high_resolution_clock::now();
            auto plain1 = encoder.encode(i);
            auto plain2 = encoder.encode(i + 1);
            auto time_encoded = chrono::high_resolution_clock::now();
            auto enc1 = encryptor.encrypt(plain1);
            auto enc2 = encryptor.encrypt(plain2);
            auto time_encrypted = chrono::high_resolution_clock::now();
            auto enc_prod = evaluator.multiply(enc1, enc2);
            auto time_multiplied = chrono::high_resolution_clock::now();
            auto enc_square = evaluator.square(enc1);
            auto time_squared = chrono::high_resolution_clock::now();
            auto enc_relin_prod = evaluator.relinearize(enc_prod);
            auto time_relinearized = chrono::high_resolution_clock::now();
            auto plain_prod = decryptor.decrypt(enc_relin_prod);
            auto time_decrypted = chrono::high_resolution_clock::now();
            int32_t result = encoder.decode_int32(plain_prod);
            auto time_decoded = chrono::high_resolution_clock::now();

            // Check the result
            int correct_result = i * (i + 1);
            if (result != correct_result)
            {
                cout << "Something went wrong (result " << result << " != " << correct_result << ")!" << endl;
            }

            if (i % 10 == 0 && i > 0)
            {
                cout << ".";
                cout.flush();
            }

            time_encode_sum += chrono::duration_cast<chrono::microseconds>(time_encoded - time_start);
            time_encrypt_sum += chrono::duration_cast<chrono::microseconds>(time_encrypted - time_encoded);
            time_multiply_sum += chrono::duration_cast<chrono::microseconds>(time_multiplied - time_encrypted);
            time_square_sum += chrono::duration_cast<chrono::microseconds>(time_squared - time_multiplied);
            time_relinearize_sum += chrono::duration_cast<chrono::microseconds>(time_relinearized - time_squared);
            time_decrypt_sum += chrono::duration_cast<chrono::microseconds>(time_decrypted - time_relinearized);
            time_decode_sum += chrono::duration_cast<chrono::microseconds>(time_decoded - time_decrypted);
        }

        cout << " done." << endl << endl;
        cout.flush();

        auto avg_encode = time_encode_sum.count() / (2 * count);
        auto avg_encrypt = time_encrypt_sum.count() / (2 * count);
        auto avg_multiply = time_multiply_sum.count() / count;
        auto avg_square = time_square_sum.count() / count;
        auto avg_relinearize = time_relinearize_sum.count() / count;
        auto avg_decrypt = time_decrypt_sum.count() / count;
        auto avg_decode = time_decode_sum.count() / count;

        cout << "Average encode: " << avg_encode << " microseconds" << endl;
        cout << "Average encrypt: " << avg_encrypt << " microseconds" << endl;
        cout << "Average multiply: " << avg_multiply << " microseconds" << endl;
        cout << "Average square: " << avg_square << " microseconds" << endl;
        cout << "Average relinearize: " << avg_relinearize << " microseconds" << endl;
        cout << "Average decrypt: " << avg_decrypt << " microseconds" << endl;
        cout << "Average decode: " << avg_decode << " microseconds" << endl;
    };

    EncryptionParameters parms;
    parms.set_poly_modulus("1x^32768 + 1");
    parms.set_coeff_modulus(ChooserEvaluator::default_parameter_options().at(32768));
    parms.set_plain_modulus(1 << 10);
    parms.set_decomposition_bit_count(64);
    SEALContext context(parms);
    performance_test(context);
    cout << endl;
}

void example_mt_timing(int th_count)
{
    print_example_banner("Example: Timing of basic operations with " + to_string(th_count) + " threads");

    auto performance_test_mt = [](int th_index, mutex &mtx, const MemoryPoolHandle &pool)
    {
        EncryptionParameters parms;
        parms.set_poly_modulus("1x^4096 + 1");
        parms.set_coeff_modulus(ChooserEvaluator::default_parameter_options().at(4096));
        parms.set_plain_modulus(1 << 10);
        parms.set_decomposition_bit_count(32);
        SEALContext context(parms, pool);

        auto poly_modulus = context.poly_modulus();
        auto coeff_modulus = context.coeff_modulus();
        auto plain_modulus = context.plain_modulus();

        KeyGenerator keygen(context, pool);
        keygen.generate(1);
        auto secret_key = keygen.secret_key();
        auto public_key = keygen.public_key();
        auto evk = keygen.evaluation_keys();

        Encryptor encryptor(context, public_key, pool);
        Decryptor decryptor(context, secret_key, pool);
        Evaluator evaluator(context, evk, pool);
        IntegerEncoder encoder(plain_modulus, 2, pool);

        chrono::microseconds time_encode_sum(0);
        chrono::microseconds time_encrypt_sum(0);
        chrono::microseconds time_multiply_sum(0);
        chrono::microseconds time_square_sum(0);
        chrono::microseconds time_relinearize_sum(0);
        chrono::microseconds time_decrypt_sum(0);
        chrono::microseconds time_decode_sum(0);

        while (!mtx.try_lock());
        cout << "Thread " << th_index << endl;
        cout << "| Memory pool address: " << hex << "0x" << &pool.operator seal::util::MemoryPool &() << dec << endl;
        cout << "| Encryption parameters:" << endl;
        cout << "| | poly_modulus: " << poly_modulus.to_string() << endl;
        cout << "| | coeff_modulus: " << coeff_modulus.to_string() << endl;
        cout << "| | plain_modulus: " << plain_modulus.to_dec_string() << endl;
        cout << "| | decomposition_bit_count: " << context.decomposition_bit_count() << endl;
        cout << "| | noise_standard_deviation: " << context.noise_standard_deviation() << endl;
        cout << "| | noise_max_deviation: " << context.noise_max_deviation() << endl;
        cout << "| Running tests ..." << endl << endl;
        cout.flush();
        mtx.unlock();

        int count = 100;

        for (int i = 0; i < count; i++)
        {
            auto time_start = chrono::high_resolution_clock::now();
            auto plain1 = encoder.encode(i);
            auto plain2 = encoder.encode(i + 1);
            auto time_encoded = chrono::high_resolution_clock::now();
            auto enc1 = encryptor.encrypt(plain1);
            auto enc2 = encryptor.encrypt(plain2);
            auto time_encrypted = chrono::high_resolution_clock::now();
            auto enc_prod = evaluator.multiply(enc1, enc2);
            auto time_multiplied = chrono::high_resolution_clock::now();
            auto enc_square = evaluator.square(enc1);
            auto time_squared = chrono::high_resolution_clock::now();
            auto enc_relin_prod = evaluator.relinearize(enc_prod);
            auto time_relinearized = chrono::high_resolution_clock::now();
            auto plain_prod = decryptor.decrypt(enc_relin_prod);
            auto time_decrypted = chrono::high_resolution_clock::now();
            int32_t result = encoder.decode_int32(plain_prod);
            auto time_decoded = chrono::high_resolution_clock::now();

            // Check the result
            int correct_result = i * (i + 1);
            if (result != correct_result)
            {
                while (!mtx.try_lock());
                cout << "Thread " << th_index << endl;
                cout << "| Something went wrong (result " << result << " != " << correct_result << ")!" << endl;
                cout.flush();
                mtx.unlock();
            }

            time_encode_sum += chrono::duration_cast<chrono::microseconds>(time_encoded - time_start);
            time_encrypt_sum += chrono::duration_cast<chrono::microseconds>(time_encrypted - time_encoded);
            time_multiply_sum += chrono::duration_cast<chrono::microseconds>(time_multiplied - time_encrypted);
            time_square_sum += chrono::duration_cast<chrono::microseconds>(time_squared - time_multiplied);
            time_relinearize_sum += chrono::duration_cast<chrono::microseconds>(time_relinearized - time_squared);
            time_decrypt_sum += chrono::duration_cast<chrono::microseconds>(time_decrypted - time_relinearized);
            time_decode_sum += chrono::duration_cast<chrono::microseconds>(time_decoded - time_decrypted);
        }

        auto avg_encode = time_encode_sum.count() / (2 * count);
        auto avg_encrypt = time_encrypt_sum.count() / (2 * count);
        auto avg_multiply = time_multiply_sum.count() / count;
        auto avg_square = time_square_sum.count() / count;
        auto avg_relinearize = time_relinearize_sum.count() / count;
        auto avg_decrypt = time_decrypt_sum.count() / count;
        auto avg_decode = time_decode_sum.count() / count;

        while (!mtx.try_lock());
        cout << "Thread " << th_index << endl;
        cout << "| Average encode: " << avg_encode << " microseconds" << endl;
        cout << "| Average encrypt: " << avg_encrypt << " microseconds" << endl;
        cout << "| Average multiply: " << avg_multiply << " microseconds" << endl;
        cout << "| Average square: " << avg_square << " microseconds" << endl;
        cout << "| Average relinearize: " << avg_relinearize << " microseconds" << endl;
        cout << "| Average decrypt: " << avg_decrypt << " microseconds" << endl;
        cout << "| Average decode: " << avg_decode << " microseconds" << endl << endl;
        cout.flush();
        mtx.unlock();
    };

    mutex print_mtx;
    vector<thread> th_vector;
    for (int i = 0; i < th_count; i++)
    {
        // Test thread-local memory pools
        //th_vector.emplace_back(performance_test_mt, i + 1, ref(print_mtx), MemoryPoolHandle::acquire_new());

        // Test global memory pool
        th_vector.emplace_back(performance_test_mt, i + 1, ref(print_mtx), MemoryPoolHandle::acquire_global());
    }
    for (int i = 0; i < th_vector.size(); i++)
    {
        th_vector[i].join();
    }
}

void example_batching_timing()
{
    print_example_banner("Example: Estimate batching performance");

    // Create encryption parameters
    EncryptionParameters parms;

    /*
    For PolyCRTBuilder it is necessary to have plain_modulus be a prime number congruent to 1 modulo
    2*degree(poly_modulus). We can use for example the following parameters:
    */
    parms.set_poly_modulus("1x^4096 + 1");
    parms.set_coeff_modulus(ChooserEvaluator::default_parameter_options().at(4096));
    parms.set_plain_modulus(40961);
    SEALContext context(parms);

    // Create the PolyCRTBuilder
    PolyCRTBuilder crtbuilder(context);
    int slot_count = crtbuilder.get_slot_count();

    cout << "Encryption parameters allow " << slot_count << " slots." << endl;

    // Create a vector of values that are to be stored in the slots. We initialize all values to 0 at this point.
    vector<BigUInt> values(slot_count, BigUInt(context.plain_modulus().bit_count(), static_cast<uint64_t>(0)));

    chrono::microseconds normal_time_compose_sum(0);
    chrono::microseconds normal_time_decompose_sum(0);
    int count = 100;
    for (int k = 0; k < count; k++)
    {
        // Set the entries to be random numbers
        random_device rd;
        for (int i = 0; i < slot_count; i++)
        {
            values[i] = (static_cast<uint64_t>(rd()) || (static_cast<uint64_t>(rd()) << 32));
        }
        auto time_start = chrono::high_resolution_clock::now();
        Plaintext plain_composed_poly = crtbuilder.compose(values);
        auto time_encoded = chrono::high_resolution_clock::now();
        crtbuilder.decompose(plain_composed_poly, values);
        auto time_decoded = chrono::high_resolution_clock::now();

        normal_time_compose_sum += chrono::duration_cast<chrono::microseconds>(time_encoded - time_start);
        normal_time_decompose_sum += chrono::duration_cast<chrono::microseconds>(time_decoded - time_encoded);
    }

    auto normal_avg_compose = normal_time_compose_sum.count() / count;
    auto normal_avg_decompose = normal_time_decompose_sum.count() / count;

    cout << "[Normal batching] Average compose: " << normal_avg_compose << " microseconds" << endl;
    cout << "[Normal batching] Average decompose: " << normal_avg_decompose << " microseconds" << endl;


    std::shared_ptr<util::ExRing> ring = util::ExRing::acquire_ring(string("A001"), 1, string("1x^1"));
    ring->init_frobe_table();
    util::ExPolyCRTBuilder exbuilder(ring, 12);

    chrono::microseconds generalized_time_compose_sum(0);
    chrono::microseconds generalized_time_decompose_sum(0);
    for (int k = 0; k < count; k++)
    {
        vector<util::ExRingElement> data(4096);
        for (int i = 0; i < 4096; i++)
            data[i] = ring->random_element();

        auto time_start = chrono::high_resolution_clock::now();
        Plaintext plain_composed_poly = exbuilder.compose(data);
        auto time_encoded = chrono::high_resolution_clock::now();
        exbuilder.decompose(plain_composed_poly, data);
        auto time_decoded = chrono::high_resolution_clock::now();

        generalized_time_compose_sum += chrono::duration_cast<chrono::microseconds>(time_encoded - time_start);
        generalized_time_decompose_sum += chrono::duration_cast<chrono::microseconds>(time_decoded - time_encoded);
    }

    auto generalized_avg_compose = generalized_time_compose_sum.count() / count;
    auto generalized_avg_decompose = generalized_time_decompose_sum.count() / count;

    cout << "[Generalized batching] Average compose: " << generalized_avg_compose << " microseconds" << endl;
    cout << "[Generalized batching] Average decompose: " << generalized_avg_decompose << " microseconds" << endl;

}


void print_example_banner(string title)
{
    if (!title.empty())
    {
        size_t title_length = title.length();
        size_t banner_length = title_length + 2 + 2 * 10;
        string banner_top(banner_length, '*');
        string banner_middle = string(10, '*') + " " + title + " " + string(10, '*');

        cout << endl
            << banner_top << endl
            << banner_middle << endl
            << banner_top << endl
            << endl;
    }
}

int64_t cpucycles()
{
    return __rdtsc();
}

void rns_variant_test()
{
    print_example_banner("Example: Full RNS Variant Scheme Test");
    chrono::microseconds time_rns_sum(0);

    RNSEncryptionParameters parms;
    parms.set_poly_modulus("1x^256 + 1");

    parms.set_plain_modulus(1 << 10);
    //this should be divide 64
    parms.set_decomposition_bit_count(8);
    vector<SmallModulus> coeff_mod_array;
	//coeff_mod_array.emplace_back(SmallModulus());
	//coeff_mod_array[0] = SmallModulus(1073479681);
    for (int i = 0; i < 2; ++i)
    {
        coeff_mod_array.emplace_back(SmallModulus());
        coeff_mod_array[i] = small_mods60[i];
    }
    parms.set_coeff_modulus(coeff_mod_array);

    RNSContext context(parms);
    IntegerEncoder encoder(context.plain_modulus(), 3);

    int input_value = 4;
    int input_value2 = 5;

    Plaintext plain_input = encoder.encode(input_value);
    Plaintext plain_input2 = encoder.encode(input_value2);
    cout << "Encoded " << input_value << " as polynomial " << plain_input.to_string() << endl;
    cout << "Encoded " << input_value2 << " as polynomial " << plain_input2.to_string() << endl;
    // Generating Keys
    cout << "Generating Keys................";
    RNSKeyGenerator generator(context);
    generator.generate(3);
    PublicKey public_key = generator.public_key();
    SecretKey secret_key = generator.secret_key();
    RNSEvaluationKeys evaluation_keys = generator.evaluation_keys();
    cout << "done." << endl;
    
    cout << "Encrypting.................";
    RNSEncryptor encryptor(context, public_key);

    Ciphertext encrypted1 = encryptor.rns_encrypt(plain_input);
    Ciphertext encrypted2 = encryptor.rns_encrypt(plain_input2);
    cout << "done." << endl;
    
    RNSEvaluator evaluator(context, evaluation_keys);
    /////////////// Homomorphic Multiplication ////////////////////
    auto time_start_mult = chrono::high_resolution_clock::now();
    Ciphertext result = evaluator.add(encrypted1, encrypted2);
    auto time_end_mult = chrono::high_resolution_clock::now();
    time_rns_sum = chrono::duration_cast<chrono::microseconds>(time_end_mult - time_start_mult);
    cout << "RNS mult takes : " << time_rns_sum.count() << " microseconds" << endl;
    
    //result = evaluator.multiply(result, result);
    //Ciphertext result = evaluator.square(encrypted1);
    //Ciphertext result = evaluator.multiply_plain(encrypted1, encoder.encode(20));
    //evaluator.transform_to_ntt(plain_input);
    //Plaintext plain_ntt = encoder.encode(20);
    //evaluator.transform_to_ntt(plain_ntt);
    //evaluator.transform_to_ntt(encrypted1);
    //evaluator.transform_from_ntt(plain_ntt);
    //Ciphertext result = evaluator.multiply_plain(encrypted1, plain_ntt);
    //evaluator.transform_from_ntt(result);
    cout << "Decrypting.................";
    RNSDecryptor decryptor(context, secret_key);
    Plaintext plain_result = decryptor.rns_decrypt(result);
    cout << "done." << endl;
    cout << "ptxt: " << plain_result.to_string() << endl;
    cout << "Result: " << encoder.decode_int64(plain_result) << endl;
    cout << "Noise: " << decryptor.invariant_noise_budget(result) << endl;
}

void fv_test_with_rns()
{
    print_example_banner("Example: FV Scheme Test");
    chrono::microseconds time_fv_sum(0);

    EncryptionParameters parms;
    parms.set_poly_modulus("1x^32768 + 1");
    parms.set_plain_modulus(1 << 10);
    parms.set_decomposition_bit_count(64);
    parms.set_coeff_modulus(ChooserEvaluator::default_parameter_options().at(32768));

    SEALContext context(parms);
    IntegerEncoder encoder(context.plain_modulus(), 3);
    int input_value = 4;
    int input_value2 = 5;
    Plaintext plain_input = encoder.encode(input_value);
    Plaintext plain_input2 = encoder.encode(input_value2);
    cout << "Encoded " << input_value << " as polynomial " << plain_input.to_string() << endl;	
    cout << "Encoded " << input_value2 << " as polynomial " << plain_input2.to_string() << endl;


    // Generating Keys
    cout << "Generating Keys............        ";
    KeyGenerator generator(context);
    generator.generate(1);
    PublicKey public_key = generator.public_key();
    SecretKey secret_key = generator.secret_key();
    cout << "done." << endl;
    
    cout << "Encrypting............        ";
    EvaluationKeys evaluation_keys = generator.evaluation_keys();
    Encryptor encryptor(context, public_key);
    Ciphertext encrypted1 = encryptor.encrypt(plain_input);
    Ciphertext encrypted2 = encryptor.encrypt(plain_input2);

    cout << "done." << endl;
    Evaluator evaluator(context, evaluation_keys);
    
    auto time_mul_start = chrono::high_resolution_clock::now();
    //Ciphertext result = evaluator.multiply_plain(encrypted1, encoder.encode(20));
    auto time_mul_end = chrono::high_resolution_clock::now();
    time_fv_sum = chrono::duration_cast<chrono::microseconds>(time_mul_end - time_mul_start);
    cout << "FV mult takes : " << time_fv_sum.count() << " microseconds" << endl;
    Ciphertext result2 = evaluator.sub(encrypted2, encrypted1);
    auto time_start = chrono::high_resolution_clock::now();
    //evaluator.relinearize(result, result, 2);
    auto time_end = chrono::high_resolution_clock::now();
    time_fv_sum = chrono::duration_cast<chrono::microseconds>(time_end - time_start);
    cout << "FV relin takes : " << time_fv_sum.count() << " microseconds" << endl;

    cout << "Decrypting.................";
    Decryptor decryptor(context, secret_key);
    Plaintext plain_result = decryptor.decrypt(result2);
    cout << "done." << endl;
    
    cout << "Result: " << encoder.decode_int64(plain_result) << endl;
    cout << "Noise: " << decryptor.invariant_noise_budget(result2) << endl;
}

void example_timing_rns()
{
	print_example_banner("Example: Timing of basic operations RNS");

	auto performance_test = [](RNSContext &context)
	{
		auto poly_modulus = context.poly_modulus();
		auto coeff_modulus = context.coeff_modulus();
		auto plain_modulus = context.plain_modulus();

		RNSKeyGenerator keygen(context);
		keygen.generate(1);
		auto secret_key = keygen.secret_key();
		auto public_key = keygen.public_key();
		auto evk = keygen.evaluation_keys();

		RNSEncryptor encryptor(context, public_key);
		RNSDecryptor decryptor(context, secret_key);
		RNSEvaluator evaluator(context, evk);
		IntegerEncoder encoder(plain_modulus);

		chrono::microseconds time_encode_sum(0);
		chrono::microseconds time_encrypt_sum(0);
		chrono::microseconds time_multiply_sum(0);
		chrono::microseconds time_square_sum(0);
		chrono::microseconds time_relinearize_sum(0);
		chrono::microseconds time_decrypt_sum(0);
		chrono::microseconds time_decode_sum(0);

		cout << "Encryption parameters:" << endl;
		cout << "| poly_modulus: " << poly_modulus.to_string() << endl;
		cout << "| coeff_modulus: " << coeff_modulus.to_string() << endl;
		cout << "| plain_modulus: " << plain_modulus.to_dec_string() << endl;
		cout << "| decomposition_bit_count: " << context.decomposition_bit_count() << endl;
		cout << "| noise_standard_deviation: " << context.noise_standard_deviation() << endl;
		cout << "| noise_max_deviation: " << context.noise_max_deviation() << endl << endl;

		int count = 10;

		cout << "Running tests ";
		for (int i = 0; i < count; i++)
		{
			auto time_start = chrono::high_resolution_clock::now();
			auto plain1 = encoder.encode(i);
			auto plain2 = encoder.encode(i + 1);
			auto time_encoded = chrono::high_resolution_clock::now();
			auto enc1 = encryptor.rns_encrypt(plain1);
			auto enc2 = encryptor.rns_encrypt(plain2);
			auto time_encrypted = chrono::high_resolution_clock::now();
			auto enc_prod = evaluator.multiply(enc1, enc2);
			auto time_multiplied = chrono::high_resolution_clock::now();
			auto enc_square = evaluator.square(enc1);
			auto time_squared = chrono::high_resolution_clock::now();
			auto enc_relin_prod = evaluator.relinearize(enc_prod);
			auto time_relinearized = chrono::high_resolution_clock::now();
			auto plain_prod = decryptor.rns_decrypt(enc_relin_prod);
			auto time_decrypted = chrono::high_resolution_clock::now();
			int32_t result = encoder.decode_int32(plain_prod);
			auto time_decoded = chrono::high_resolution_clock::now();

			// Check the result
			int correct_result = i * (i + 1);
			if (result != correct_result)
			{
				cout << "Something went wrong (result " << result << " != " << correct_result << ")!" << endl;
			}

			if (i % 10 == 0 && i > 0)
			{
				cout << ".";
				cout.flush();
			}

			time_encode_sum += chrono::duration_cast<chrono::microseconds>(time_encoded - time_start);
			time_encrypt_sum += chrono::duration_cast<chrono::microseconds>(time_encrypted - time_encoded);
			time_multiply_sum += chrono::duration_cast<chrono::microseconds>(time_multiplied - time_encrypted);
			time_square_sum += chrono::duration_cast<chrono::microseconds>(time_squared - time_multiplied);
			time_relinearize_sum += chrono::duration_cast<chrono::microseconds>(time_relinearized - time_squared);
			time_decrypt_sum += chrono::duration_cast<chrono::microseconds>(time_decrypted - time_relinearized);
			time_decode_sum += chrono::duration_cast<chrono::microseconds>(time_decoded - time_decrypted);
		}

		cout << " done." << endl << endl;
		cout.flush();

		auto avg_encode = time_encode_sum.count() / (2 * count);
		auto avg_encrypt = time_encrypt_sum.count() / (2 * count);
		auto avg_multiply = time_multiply_sum.count() / count;
		auto avg_square = time_square_sum.count() / count;
		auto avg_relinearize = time_relinearize_sum.count() / count;
		auto avg_decrypt = time_decrypt_sum.count() / count;
		auto avg_decode = time_decode_sum.count() / count;

		cout << "Average encode: " << avg_encode << " microseconds" << endl;
		cout << "Average encrypt: " << avg_encrypt << " microseconds" << endl;
		cout << "Average multiply: " << avg_multiply << " microseconds" << endl;
		cout << "Average square: " << avg_square << " microseconds" << endl;
		cout << "Average relinearize: " << avg_relinearize << " microseconds" << endl;
		cout << "Average decrypt: " << avg_decrypt << " microseconds" << endl;
		cout << "Average decode: " << avg_decode << " microseconds" << endl;
	};

	RNSEncryptionParameters parms;
	parms.set_poly_modulus("1x^16384 + 1");
	vector<SmallModulus> coeff_mod_array;
	for (int i = 0; i < 14; ++i)
	{
		coeff_mod_array.emplace_back(SmallModulus());
		coeff_mod_array[i] = small_mods60[i];
	}
	
	/*coeff_mod_array.emplace_back(1073479681);
	coeff_mod_array.emplace_back(1072496641);
	coeff_mod_array.emplace_back(1071513601);
	*///coeff_mod_array[14] = SmallModulus(1073479681);
	parms.set_coeff_modulus(coeff_mod_array);
	parms.set_plain_modulus(1 << 10);
	parms.set_decomposition_bit_count(60);
	RNSContext context(parms);
	performance_test(context);
	cout << endl;
}

void test_for_added_functions()
{
	RNSEncryptionParameters parms;
	parms.set_poly_modulus("1x^32 + 1");
	vector<SmallModulus> coeff_mod_array;
	for (int i = 0; i < 12; ++i)
	{
		coeff_mod_array.emplace_back(SmallModulus());
		coeff_mod_array[i] = small_mods60[i];
	}
	parms.set_coeff_modulus(coeff_mod_array);
	parms.set_plain_modulus(24137569);
	parms.set_plain_modulus_base(17);
	parms.set_plain_modulus_exponent(6);
	parms.set_decomposition_bit_count(60);
	RNSContext context(parms);

	auto poly_modulus = context.poly_modulus();
	auto coeff_modulus = context.coeff_modulus();
	auto plain_modulus = context.plain_modulus();

	RNSKeyGenerator keygen(context);
	keygen.generate(1, 64, true);

	string temp_string = "1x^" + to_string(321);
	Plaintext power_of_X = temp_string;
	cout << power_of_X.to_string() << endl;

	auto secret_key = keygen.secret_key();
	auto public_key = keygen.public_key();
	auto evk = keygen.evaluation_keys();
	
	RNSEncryptor encryptor(context, public_key);
	RNSDecryptor decryptor(context, secret_key);
	RNSEvaluator evaluator(context, evk);
	evaluator.set_galois_keys(keygen.galois_keys());

	Plaintext plain("3x^16 + 12x^4 + 11x^3 + 1x^2 + 1x^1 + 123");
	Plaintext change("1x^31");
	{
		Plaintext const_plain("164cce6"); // inverse_mod(n, plain_mod)
		int n = parms.poly_modulus().coeff_count() - 1;
		int m = n << 1;
		vector<uint64_t> galois_elts;
		for (int i = 0; i < 5; i++)
		{
			galois_elts.push_back((n + power(2, i)) / power(2, i));
		}

		Ciphertext cipher = encryptor.rns_encrypt(plain);
		cipher = evaluator.multiply_plain(cipher, change);
		cipher = evaluator.negate(cipher);
		Ciphertext result = cipher;
		
		for (int i = 0; i < 5; i++)
		{
			Ciphertext temp;
			temp = evaluator.apply_galois(result, galois_elts[i]);
			result = evaluator.add(temp, result);
		}
		result = evaluator.multiply_plain(result, const_plain);
		cout << decryptor.rns_decrypt(result).to_string() << endl;
	}
}

void test_for_recryptor()
{
	RNSEncryptionParameters parms;
	parms.set_poly_modulus("1x^2048 + 1");
	vector<SmallModulus> coeff_mod_array;
	int bit_count = 0;
	for (int i = 0; i < 8; ++i)
	{
		coeff_mod_array.emplace_back(SmallModulus());
		coeff_mod_array[i] = small_mods60[i];
		bit_count += coeff_mod_array[i].bit_count();
	}
	cout << "Total bit count = " << bit_count << endl;
	
	parms.set_coeff_modulus(coeff_mod_array);
	parms.set_plain_modulus(257);
	parms.set_plain_modulus_base(257);
	parms.set_plain_modulus_exponent(1);
	parms.set_decomposition_bit_count(60);
	RNSContext context(parms);

	auto poly_modulus = context.poly_modulus();
	auto coeff_modulus = context.coeff_modulus();
	auto plain_modulus = context.plain_modulus();

	RNSKeyGenerator keygen(context);
	keygen.generate(1, 64, true);

	auto secret_key = keygen.secret_key();
	auto public_key = keygen.public_key();
	auto evk = keygen.evaluation_keys();
	
	RNSEncryptor encryptor(context, public_key);
	RNSDecryptor decryptor(context, secret_key);
	RNSEvaluator evaluator(context, evk);
 
  evaluator.set_galois_keys(keygen.galois_keys());
	
	auto time_recryptor_start = chrono::high_resolution_clock::now();

	//RNSRecryptor recryptor(context, secret_key, evk, "1x^32 + 1FCx^16 + 2D8", 6, true); // n = 64 p = 3 r = 2
	//RNSRecryptor recryptor(context, secret_key, evk, "1x^128 + 1FCx^64 + 2D8", 6, true); // n = 64 p = 3 r = 2

	//RNSRecryptor recryptor(context, secret_key, evk, "1x^2 + 2C3x^1 + 3F00", 2); // n = 128 p = 127 r = 1


	// n = 128, p = 31, e = 3
	//Recryptor recryptor(context, secret_key, evaluation_key, "1x^8 + 1699x^4 + 745E", 3, new_coeff_modulus_);		


	// n = 128, p = 127, e = 2
	//RNSRecryptor recryptor(context, secret_key, evk, "1x^2 + 2C3x^1 + 3F00", 2);

	// n = 1024, p = 127, e = 2
	//RNSRecryptor recryptor(context, secret_key, evk, "1x^16 + 2C3x^8 + 3F00",2); // n = 1024, p = 127, e = 2

	// n = 2048, p = 127, e = 2																															   
	//RNSRecryptor recryptor(context, secret_key, evk, "1x^32 + 2C3x^16 + 3F00", 2);

	// n = 2048, p = 257, e = 2																															   
	RNSRecryptor recryptor(context, secret_key, evk, keygen.galois_keys(), "1x^16 + 869A", 2);

	// n = 4096, p = 257, e = 2																															   
	//RNSRecryptor recryptor(context, secret_key, evk, "1x^32 + 869A", 2);

	// n = 4096, p = 383, e = 2																															   
	//Recryptor recryptor(context, secret_key, evaluation_key, "1x^64 + 17256x^32 + 23D00", 2, new_coeff_modulus_);

	// n = 4096, p = 769, e = 2
	//RNSRecryptor recryptor(context, secret_key, evk, "1x^32 + 7A2BD", 2);



	// n = 8192, p = 257, e = 2
	//RNSRecryptor recryptor(context, secret_key, evk, "1x^64 + 869A", 2);

	// n = 16384, p = 257, e = 2
	//RNSRecryptor recryptor(context, secret_key, evk, "1x^128 + 869A", 2);


	// n = 16384, p = 127, e = 2
	//Recryptor recryptor(context, secret_key, evaluation_key, "1x^256 + 2C3x^128 + 3F00", 2, new_coeff_modulus_);

	// n = 16384, p = 383, e = 2
	//Recryptor recryptor(context, secret_key, evaluation_key, "1x^256 + 17256x^128 + 23D00", 2, new_coeff_modulus_);


	// n = 16384, p = 769, e = 2
	//Recryptor recryptor(context, secret_key, evaluation_key, "1x^128 + 7A2BD", 2, new_coeff_modulus_);

	// RNSRecryptor recryptor(context, secret_key, evk, "1x^16 + 2C3x^8 + 3F00", 2); // n = 1024 p = 127 r = 1

	// n = 16384, p = 257, e = 2
	//RNSRecryptor recryptor(context, secret_key, evk, "1x^128 + 869A", 2);


	auto time_recryptor_end = chrono::high_resolution_clock::now();
	cout << "recrypt initializing time : " << (double)(chrono::duration_cast<chrono::microseconds>(time_recryptor_end - time_recryptor_start).count()) / 1000000. << " seconds" << endl;

	Plaintext plain("1x^2 + 1x^1 + 1");
  auto time_encrypt_start = chrono::high_resolution_clock::now();
	Ciphertext cipher = encryptor.rns_encrypt(plain);
	auto time_recrypt_start = chrono::high_resolution_clock::now();
	Ciphertext cipher_new = recryptor.recrypt(cipher);
	auto time_recrypt_end = chrono::high_resolution_clock::now();

	cout << "Before recrypt : " << decryptor.rns_decrypt(cipher).to_string() << endl;
	cout << "After recrypt : " << decryptor.rns_decrypt(cipher_new).to_string() << endl;
	cout << "After recrypt noise : " << decryptor.invariant_noise_budget(cipher_new) << endl;
  cout << "encrypt time : " << (double)(chrono::duration_cast<chrono::microseconds>(time_recrypt_start - time_encrypt_start).count()) / 1000000. << " seconds" << endl;
	cout << "recrypt time : " << (double)(chrono::duration_cast<chrono::microseconds>(time_recrypt_end - time_recrypt_start).count()) / 1000000. << " seconds" << endl;

}

void example_floor()
{
	RNSEncryptionParameters parms;
	parms.set_poly_modulus("1x^256 + 1");
	vector<SmallModulus> coeff_mod_array;
	for (int i = 0; i < 16; ++i)
	{
		coeff_mod_array.emplace_back(SmallModulus());
		coeff_mod_array[i] = small_mods60[i];
	}
	parms.set_coeff_modulus(coeff_mod_array);

	// plaintext modu = p^3. 
	parms.set_plain_modulus_base(127);
	parms.set_plain_modulus(2048383);
	parms.set_plain_modulus_exponent(3);
	parms.set_decomposition_bit_count(60);
	RNSContext context(parms);

	auto poly_modulus = context.poly_modulus();
	auto coeff_modulus = context.coeff_modulus();
	auto plain_modulus = context.plain_modulus();

	cout << "Keygen....";
	RNSKeyGenerator keygen(context);
	keygen.generate(1, 20, true);
	cout << "Done!" << endl;

	auto secret_key = keygen.secret_key();
	auto public_key = keygen.public_key();
	auto evk = keygen.evaluation_keys();

	RNSEncryptor encryptor(context, public_key);
	RNSDecryptor decryptor(context, secret_key);
	RNSEvaluator evaluator(context, evk);
	//evaluator.set_galois_keys(keygen.galois_keys());

	cout << "Recrypt related keygen...";
	RNSRecryptor recryptor(context, secret_key, evk,  keygen.galois_keys(), "1x^4 + 5607f85x^2 + f817e00", 4, true); //p=127 e=4 n=256
	//RNSRecryptor recryptor(context, secret_key, evk, "1x^8 + 5607f85x^4 + f817e00", 4, true); //p=127 e=4 n=512
	//RNSRecryptor recryptor(context, secret_key, evk, keygen.galois_keys(), "1x^16 + 5607f85x^8 + f817e00", 4, true); //p=127 e=4 n=1024
	//RNSRecryptor recryptor(context, secret_key, evk, "1x^32 + 5607f85x^16 + f817e00", 4, true); //p=127 e=4 n=2048
	//RNSRecryptor recryptor(context, secret_key, evk, "1x^64 + 5607f85x^32 + f817e00", 4, true); //p=127 e=4 n=4096
	//RNSRecryptor recryptor(context, secret_key, evk, "1x^128 + 5607f85x^64 + f817e00", 4, true); //p=127 e=4 n=8192
	//RNSRecryptor recryptor(context, secret_key, evk, "1x^512 + 5607f85x^256 + f817e00", 4, true); //p=127 e=4 n=32768
	cout << "done" << endl;

	std::shared_ptr<ExRing> ring = ExRing::acquire_ring(string("7f"), 3, string("1x^4 + 179afbx^2 + 1f417e")); // p=127 e=3 n=256
	//std::shared_ptr<ExRing> ring = ExRing::acquire_ring(string("7f"), 3, string("1x^8 + 179afbx^4 + 1f417e")); // p=127 e=3 n=512
	//std::shared_ptr<ExRing> ring = ExRing::acquire_ring(string("7f"), 3, string("1x^16 + 179afbx^8 + 1f417e")); // p=127 e=3 n=1024
	//std::shared_ptr<ExRing> ring = ExRing::acquire_ring(string("7f"), 3, string("1x^32 + 179afbx^16 + 1f417e")); // p=127 e=3 n=2048
	//std::shared_ptr<ExRing> ring = ExRing::acquire_ring(string("7f"), 3, string("1x^64 + 179afbx^32 + 1f417e")); // p=127 e=3 n=4096
	//std::shared_ptr<ExRing> ring = ExRing::acquire_ring(string("7f"), 3, string("1x^128 + 179afbx^64 + 1f417e")); // p=127 e=3 n=8192
	//std::shared_ptr<ExRing> ring = ExRing::acquire_ring(string("7f"), 3, string("1x^512 + 179afbx^256 + 1f417e")); // p=127 e=3 n=32768

	int m = (parms.poly_modulus().coeff_count() - 1) << 1;
	ring->init_frobe_table(m);
	ExPolyCRTBuilder builder(ring, get_power_of_two(parms.poly_modulus().coeff_count() - 1));
	int numslot = builder.slot_count();
	vector<ExRingElement> vec_exringelt(numslot);
	BigUInt input("9F"); // 127 + 32. 
	for (int i = 0; i < numslot; i++)
	{
		vec_exringelt[i] = ExRingElement(ring, input.to_string());
		// input += 127; 
	}
	Plaintext plain = builder.compose(vec_exringelt);
	Ciphertext cipher = encryptor.rns_encrypt(plain);
	Ciphertext result = recryptor.floor_simd(cipher, 1);

	cout << "Msg before floor function: ";
	for (int i = 0; i < numslot; i++)
	{
		cout << vec_exringelt[i].to_string() << ", ";
	}

	cout << "Msg after floor function: ";
	for (int i = 0; i < numslot; i++)
	{
		cout <<  builder.decompose(decryptor.rns_decrypt(result))[i].to_string() << ", ";
	}
	cout << endl;
	cout << "Noise after floor function: " << decryptor.invariant_noise_budget(result) << endl;
}


void example_sigmoid()
{
	RNSEncryptionParameters parms;
	parms.set_poly_modulus("1x^256 + 1");
	vector<SmallModulus> coeff_mod_array;
	for (int i = 0; i < 16; ++i)
	{
		coeff_mod_array.emplace_back(SmallModulus());
		coeff_mod_array[i] = small_mods60[i];
	}
	parms.set_coeff_modulus(coeff_mod_array);

	int p = 127;
	int e = 3;
	uint64_t ptoe = power(p, e);

	// plaintext modu = p^3. 
	parms.set_plain_modulus_base(127);
	parms.set_plain_modulus(ptoe);
	parms.set_plain_modulus_exponent(e);
	parms.set_decomposition_bit_count(60);
	RNSContext context(parms);

	auto poly_modulus = context.poly_modulus();
	auto coeff_modulus = context.coeff_modulus();
	auto plain_modulus = context.plain_modulus();

	cout << "Keygen....";
	RNSKeyGenerator keygen(context);
	keygen.generate(1, 64, true);
	cout << "Done!" << endl;

	auto secret_key = keygen.secret_key();
	auto public_key = keygen.public_key();
	auto evk = keygen.evaluation_keys();

	RNSEncryptor encryptor(context, public_key);
	RNSDecryptor decryptor(context, secret_key);
	RNSEvaluator evaluator(context, evk);
	//evaluator.set_galois_keys(keygen.galois_keys());

	cout << "Recrypt related keygen...";
	RNSRecryptor recryptor(context, secret_key, evk, keygen.galois_keys(), "1x^4 + 5607f85x^2 + f817e00" , 4, true); //p=127 e=4 n=256
	//RNSRecryptor recryptor(context, secret_key, evk, "1x^8 + 5607f85x^4 + f817e00", 4, true); //p=127 e=4 n=512
	//RNSRecryptor recryptor(context, secret_key, evk, keygen.galois_keys(), "1x^16 + 5607f85x^8 + f817e00", 4, true); //p=127 e=4 n=1024
																													 //RNSRecryptor recryptor(context, secret_key, evk, "1x^32 + 5607f85x^16 + f817e00", 4, true); //p=127 e=4 n=2048
																													 //RNSRecryptor recryptor(context, secret_key, evk, "1x^64 + 5607f85x^32 + f817e00", 4, true); //p=127 e=4 n=4096
																													 //RNSRecryptor recryptor(context, secret_key, evk, "1x^128 + 5607f85x^64 + f817e00", 4, true); //p=127 e=4 n=8192
																													 //RNSRecryptor recryptor(context, secret_key, evk, "1x^512 + 5607f85x^256 + f817e00", 4, true); //p=127 e=4 n=32768
	cout << "done" << endl;

	std::shared_ptr<ExRing> ring = ExRing::acquire_ring(string("7f"), 3, string("1x^4 + 179afbx^2 + 1f417e")); // p=127 e=3 n=256
	//std::shared_ptr<ExRing> ring = ExRing::acquire_ring(string("7f"), 3, string("1x^8 + 179afbx^4 + 1f417e")); // p=127 e=3 n=512
	//std::shared_ptr<ExRing> ring = ExRing::acquire_ring(string("7f"), 3, string("1x^16 + 179afbx^8 + 1f417e")); // p=127 e=3 n=1024


	int m = (parms.poly_modulus().coeff_count() - 1) << 1;
	ring->init_frobe_table(m);
	ExPolyCRTBuilder builder(ring, get_power_of_two(parms.poly_modulus().coeff_count() - 1));
	int numslot = builder.slot_count();
	vector<ExRingElement> vec_exringelt(numslot);
	//BigUInt input = "1F417E"; // p^3 - 1;
	BigUInt input("0");
	// vector<string> inputs = 
	for (int i = 0; i < numslot; i++)
	{
		vec_exringelt[i] = ExRingElement(ring, input.to_string());
	}
	Plaintext plain = builder.compose(vec_exringelt);
	Ciphertext cipher = encryptor.rns_encrypt(plain);
	Ciphertext result; 

	auto plain_vec = builder.decompose(decryptor.rns_decrypt(cipher));
	cout << "plain = ";
	for (int i = 0; i < numslot; i++)
	{
		cout << plain_vec[i].pointer()[0] << ", ";
	}
	cout << endl;

	recryptor.evaluate_twotimes_sigmoid(cipher, result, decryptor);

	cout << "Msg after evaluating  2*sigmoid function: ";
	//int x = -1; 
	int x = 0;
	auto result_vec = builder.decompose(decryptor.rns_decrypt(result));
	for (int i = 0; i < numslot; i++)
	{
		cout << "x = " << (double)(x) / p << endl;
		cout << "2 sigma(x) =  ";  
		int64_t value = result_vec[i].pointer()[0];
		if (value > ptoe / 2) value -= ptoe; 
		cout << (double)(value) / p << endl;
		// x -= 50;
	}
	cout << endl;
	cout << "Noise after sigmoid function: " << decryptor.invariant_noise_budget(result) << endl;
}


void example_applygalois()
{
	RNSEncryptionParameters parms;
	parms.set_poly_modulus("1x^8192 + 1");
	vector<SmallModulus> coeff_mod_array;
	for (int i = 0; i < 4; ++i)
	{
		coeff_mod_array.emplace_back(SmallModulus());
		coeff_mod_array[i] = small_mods60[i];
	}
	parms.set_coeff_modulus(coeff_mod_array);
	parms.set_plain_modulus(1 << 10);
	parms.set_decomposition_bit_count(60);
	RNSContext context(parms);

	auto poly_modulus = context.poly_modulus();
	auto coeff_modulus = context.coeff_modulus();
	auto plain_modulus = context.plain_modulus();

	cout << "Generating keys ..." << endl;
	RNSKeyGenerator keygen(context);
	keygen.generate(1, 64, true); // generate key with SKhammingweight = 64 and logn number of galois key
	cout << "... key generation complete" << endl;

	auto secret_key = keygen.secret_key();
	auto public_key = keygen.public_key();
	auto evk = keygen.evaluation_keys();

	RNSEncryptor encryptor(context, public_key);
	RNSDecryptor decryptor(context, secret_key);
	RNSEvaluator evaluator(context, evk);
	evaluator.set_galois_keys(keygen.galois_keys());

	int count = 0;
	double time_average = 0.0;
    random_device rd;
	while (count < 20)
	{
		int n = (parms.poly_modulus().coeff_count() - 1);
		uint64_t galois_elt = (2 * rd() + 1) % n;
		Plaintext plain("1x^1 + 1");
		Ciphertext cipher = encryptor.rns_encrypt(plain);
		auto time_start = chrono::high_resolution_clock::now();
		Ciphertext cipher_applygalois = evaluator.apply_galois(cipher, galois_elt);
		auto time_end = chrono::high_resolution_clock::now();
		string result = "1x^" + to_string(galois_elt) + " + 1";
		if (decryptor.rns_decrypt(cipher_applygalois).to_string() != result)
		{
			cout << "Apply Galois is not Correct" << endl;
		}
		else
		{
			time_average += (double)chrono::duration_cast<chrono::microseconds>(time_end - time_start).count() / 1000.;
		}
		count++;
	}
	cout << "Galois Average Time : " << time_average / 20 << " milliseconds" << endl;
}


void test_digitextraction() {
	RNSEncryptionParameters parms;
	parms.set_poly_modulus("1x^256 + 1");
	vector<SmallModulus> coeff_mod_array;
	int bit_count = 0;
	for (int i = 0; i < 16; ++i)
	{
		coeff_mod_array.emplace_back(SmallModulus());
		coeff_mod_array[i] = small_mods60[i];
		bit_count += coeff_mod_array[i].bit_count();
	}
	cout << "Total bit count = " << bit_count << endl;

	int p = 127; 
	int e = 4; 


	parms.set_coeff_modulus(coeff_mod_array);
	parms.set_plain_modulus(power(p,e));
	parms.set_plain_modulus_base(p);
	parms.set_plain_modulus_exponent(e);
	parms.set_decomposition_bit_count(60);
	RNSContext context(parms);

	auto poly_modulus = context.poly_modulus();
	auto coeff_modulus = context.coeff_modulus();
	auto plain_modulus = context.plain_modulus();

	RNSKeyGenerator keygen(context);
	keygen.generate(1);

	auto secret_key = keygen.secret_key();
	auto public_key = keygen.public_key();
	auto evk = keygen.evaluation_keys();

	RNSEncryptor encryptor(context, public_key);
	RNSDecryptor decryptor(context, secret_key);
	RNSEvaluator evaluator(context, evk);


	Plaintext plain("F23B9C7"); // [-60,1,0,-3]
	Ciphertext cipher = encryptor.rns_encrypt(plain);
	auto time_recrypt_start = chrono::high_resolution_clock::now();
	vector<Ciphertext> cipher_new = evaluator.DigitExtraction(cipher, p, e, e-1,  true);
	auto time_recrypt_end = chrono::high_resolution_clock::now();

	cout << "Msg before extraction: " << decryptor.rns_decrypt(cipher).to_string() << endl;
	cout << "Msg after extraction : ";
	for (int i = 0; i < cipher_new.size(); i++) {
		int64_t value = static_cast<int64_t> (decryptor.rns_decrypt(cipher_new[i])[0].pointer()[0]);
		if (value > power(p, e) / 2) value -= power(p, e); 
		// cout << "value = " << value << ", ";
		cout << i << "-th digit = " << (double) value / power(p, i);
		cout << ", noise : " << decryptor.invariant_noise_budget(cipher_new[i]) << endl;
	}

	Ciphertext test;
	test = cipher_new[1];
	evaluator.add(test, cipher_new[2]); 
	cout << " middle digits = " << decryptor.rns_decrypt(test)[0].pointer()[0] << endl;

}

void test_mulbylearningrate() {
	RNSEncryptionParameters parms;
	parms.set_poly_modulus("1x^256 + 1");
	vector<SmallModulus> coeff_mod_array;
	for (int i = 0; i < 16; ++i)
	{
		coeff_mod_array.emplace_back(SmallModulus());
		coeff_mod_array[i] = small_mods60[i];
	}
	parms.set_coeff_modulus(coeff_mod_array);

	int p = 127;
	int e = 3;
	uint64_t ptoe = power(p, e);

	// plaintext modu = p^3. 
	parms.set_plain_modulus_base(127);
	parms.set_plain_modulus(ptoe);
	parms.set_plain_modulus_exponent(e);
	parms.set_decomposition_bit_count(60);
	RNSContext context(parms);

	auto poly_modulus = context.poly_modulus();
	auto coeff_modulus = context.coeff_modulus();
	auto plain_modulus = context.plain_modulus();

	cout << "Keygen....";
	RNSKeyGenerator keygen(context);
	keygen.generate(1, 64, true);
	cout << "Done!" << endl;

	auto secret_key = keygen.secret_key();
	auto public_key = keygen.public_key();
	auto evk = keygen.evaluation_keys();

	RNSEncryptor encryptor(context, public_key);
	RNSDecryptor decryptor(context, secret_key);
	RNSEvaluator evaluator(context, evk);
	//evaluator.set_galois_keys(keygen.galois_keys());

	cout << "Recrypt related keygen...";
	RNSRecryptor recryptor(context, secret_key, evk, keygen.galois_keys(), "1x^4 + 5607f85x^2 + f817e00", 4, true); //p=127 e=4 n=256
																													//RNSRecryptor recryptor(context, secret_key, evk, "1x^8 + 5607f85x^4 + f817e00", 4, true); //p=127 e=4 n=512
																													//RNSRecryptor recryptor(context, secret_key, evk, keygen.galois_keys(), "1x^16 + 5607f85x^8 + f817e00", 4, true); //p=127 e=4 n=1024
																													//RNSRecryptor recryptor(context, secret_key, evk, "1x^32 + 5607f85x^16 + f817e00", 4, true); //p=127 e=4 n=2048
																													//RNSRecryptor recryptor(context, secret_key, evk, "1x^64 + 5607f85x^32 + f817e00", 4, true); //p=127 e=4 n=4096
																													//RNSRecryptor recryptor(context, secret_key, evk, "1x^128 + 5607f85x^64 + f817e00", 4, true); //p=127 e=4 n=8192																													//RNSRecryptor recryptor(context, secret_key, evk, "1x^512 + 5607f85x^256 + f817e00", 4, true); //p=127 e=4 n=32768
	cout << "done" << endl;

	std::shared_ptr<ExRing> ring = ExRing::acquire_ring(string("7f"), 3, string("1x^4 + 179afbx^2 + 1f417e")); // p=127 e=3 n=256
																											   //std::shared_ptr<ExRing> ring = ExRing::acquire_ring(string("7f"), 3, string("1x^8 + 179afbx^4 + 1f417e")); // p=127 e=3 n=512
																											   //std::shared_ptr<ExRing> ring = ExRing::acquire_ring(string("7f"), 3, string("1x^16 + 179afbx^8 + 1f417e")); // p=127 e=3 n=1024


	int m = (parms.poly_modulus().coeff_count() - 1) << 1;
	ring->init_frobe_table(m);
	ExPolyCRTBuilder builder(ring, get_power_of_two(parms.poly_modulus().coeff_count() - 1));
	int numslot = builder.slot_count();
	vector<ExRingElement> vec_exringelt(numslot);
	BigUInt input("1F417E"); // p^3 - 1;
	//BigUInt input = "0";
	for (int i = 0; i < numslot; i++)
	{
		vec_exringelt[i] = ExRingElement(ring, input.to_string());
		input -= 50;
	}
	Plaintext plain = builder.compose(vec_exringelt);
	Ciphertext cipher = encryptor.rns_encrypt(plain);
	Ciphertext result;

	auto plain_vec = builder.decompose(decryptor.rns_decrypt(cipher));
	cout << "plain = ";
	for (int i = 0; i < numslot; i++)
	{
		cout << plain_vec[i].pointer()[0] << ", ";
	}
	cout << endl;

	recryptor.multiply_by_learningrate(cipher, result);

	cout << "Msg after multiplying by 0.0016 : ";
	int x = -1;
	auto result_vec = builder.decompose(decryptor.rns_decrypt(result));
	for (int i = 0; i < numslot; i++)
	{
		cout << "x = " << (double)(x) / p  << ", 0.0016 x (expected) = " << (double) (x*0.0016) /p << endl;
		cout << "0.0016 x =  ";
		int64_t value = result_vec[i].pointer()[0];
		if (value > ptoe / 2) value -= ptoe;
		cout << (double)(value) / p << endl;
		x -= 50;
	}
	cout << endl;
	cout << "Noise after evaluation: " << decryptor.invariant_noise_budget(result) << endl;
}



void test_combineexpand() {
	RNSEncryptionParameters parms;
	parms.set_poly_modulus("1x^256 + 1");
	vector<SmallModulus> coeff_mod_array;
	for (int i = 0; i < 16; ++i)
	{
		coeff_mod_array.emplace_back(SmallModulus());
		coeff_mod_array[i] = small_mods60[i];
	}
	parms.set_coeff_modulus(coeff_mod_array);

	// plaintext modu = p^3. 
	parms.set_plain_modulus_base(127);
	parms.set_plain_modulus(2048383);
	parms.set_plain_modulus_exponent(3);
	parms.set_decomposition_bit_count(60);
	RNSContext context(parms);

	auto poly_modulus = context.poly_modulus();
	auto coeff_modulus = context.coeff_modulus();
	auto plain_modulus = context.plain_modulus();

	cout << "Keygen....";
	RNSKeyGenerator keygen(context);
	keygen.generate(1, 64, true);
	cout << "Done!" << endl;

	auto secret_key = keygen.secret_key();
	auto public_key = keygen.public_key();
	auto evk = keygen.evaluation_keys();

	RNSEncryptor encryptor(context, public_key);
	RNSDecryptor decryptor(context, secret_key);
	RNSEvaluator evaluator(context, evk);
	evaluator.set_galois_keys(keygen.galois_keys());

	cout << "Recrypt related keygen...";
    RNSRecryptor recryptor(context, secret_key, evk, keygen.galois_keys(), "1x^4 + 5607f85x^2 + f817e00", 4, true); //p=127 e=4 n=256
	// RNSRecryptor recryptor(context, secret_key, evk, keygen.galois_keys(), "1x^4 + 179afbx^2 + 1f417e", 4, true); //p=127 e=4 n=256
	cout << "done" << endl;

	std::shared_ptr<ExRing> ring = ExRing::acquire_ring(string("7f"), 3, string("1x^4 + 179afbx^2 + 1f417e")); // p=127 e=3 n=256

	int m = (parms.poly_modulus().coeff_count() - 1) << 1;
	ring->init_frobe_table(m);
	ExPolyCRTBuilder builder(ring, get_power_of_two(parms.poly_modulus().coeff_count() - 1));
	int numslot = builder.slot_count();
	vector<ExRingElement> vec_exringelt(numslot);
	BigUInt input("0"); 
	for (int i = 0; i < numslot; i++)
	{
		vec_exringelt[i] = ExRingElement(ring, input.to_string());
		input += 127; 
	}
	Plaintext plain = builder.compose(vec_exringelt);
	Ciphertext cipher = encryptor.rns_encrypt(plain);


	vector<Plaintext> unit_vecs;
	unit_vecs.resize(numslot);
	for (int i = 0; i < numslot; i++)
	{
		vector<ExRingElement> unit_vector(numslot);
		for (int j = 0; j < numslot; j++)
		{
			if (i == j)
			{
				unit_vector[j] = ExRingElement(ring, "1");
			}
			else
			{
				unit_vector[j] = ExRingElement(ring, "0");
			}
		}
		unit_vecs[i] = builder.compose(unit_vector);
	}


	vector<int64_t> gcd = xgcd(3, m);
	if (gcd[0] != 1)
	{
		throw "Inverse does not exist!";
	}
	if (gcd[1] < 0)
	{
		gcd[1] += m;
	}
	uint64_t three_inv_mod_m = gcd[1];

	int log_slot_num = get_power_of_two(numslot);
	vector<uint64_t> power_of_three_coset_index(log_slot_num);
	for (int i = 0; i < log_slot_num; i++)
	{
		power_of_three_coset_index[log_slot_num - i - 1] = three_inv_mod_m;
		three_inv_mod_m *= three_inv_mod_m;
		three_inv_mod_m %= m;
	}

	cout << "index done " << endl;
	
	vector<Ciphertext> result(10);

	
	evaluator.expand(cipher, result, unit_vecs, power_of_three_coset_index);

	cout << "Msg before expand function: ";
	for (int i = 0; i < numslot; i++)
	{
		cout << vec_exringelt[i].to_string() << ", ";
	}


	cout << "Msg after expand function: ";
	for (int i = 0; i < result.size(); i++)
	{
		cout << "Noise after expand:  " << decryptor.invariant_noise_budget(result[i]) << endl;
		vector<ExRingElement> plain_vec = builder.decompose(decryptor.rns_decrypt(result[i]));
		for (int j = 0; j < numslot; j++) {
			cout << plain_vec[j].to_string() << ", ";
		}
		cout << endl;
	}

	cout << endl; 

	Ciphertext combined; 
	evaluator.combine(result, combined, unit_vecs);

	cout << "Noise after combine:  " << decryptor.invariant_noise_budget(combined) << endl;
	vector<ExRingElement> plain_vec = builder.decompose(decryptor.rns_decrypt(combined));
	for (int j = 0; j < numslot; j++) {
		cout << plain_vec[j].to_string() << ", ";
	}
	cout << endl;
}


void test_centereddigits() {
	int p = 127; 
	int64_t x = p; 
	int e = 4;
	vector<int> digits = centered_digits(x, p, e);
	cout << "["; 
	for (int i = 0; i < digits.size(); i++) {
		cout << digits[i] << ", ";
	}
	cout << "]" << endl;
}
