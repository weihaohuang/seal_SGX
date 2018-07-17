#include "SealEnclaveTest_u.h"

#include "sgx_urts.h"
#include <stdio.h>
// #include <tchar.h> // windows environment
#include <string.h>
#include "Matrix.h"
#include <vector>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <sstream>
#include <chrono>
#include <unordered_map>
#include "../Seal_OutEnclave/seal/seal.h"
#include "TestData.h"
//#include "ReadData.h"
#include "MakeConfigure.h"

using namespace std;
using namespace seal;

#define ENCLAVE_FILE "SealEnclaveTest.signed.so"
#define MAX_BUF_LEN 600000
#define round(x) (x - floor(x) >= 0.5 ? floor(x) + 1 : floor(x)) 
#define SIGNIFICANT_FIGURES 2
#define Simplify_NewSigmoid
#define DEBUG
#define SHOWRESULT

//vector<vector<double>> X;
//vector<double> Y;
sgx_enclave_id_t      eid;



PublicKey          public_key;
#if defined (DEBUG) || defined(SHOWRESULT)
SecretKey               secret_key;
#endif // DEBUG

//------------------------------------------------------------------------------
void ocall_print(const char* str) {
	printf("\033[96m%s\033[0m", str);
}


//--------------< Struct for the project >---------------------------------------------
static struct Configure
{
	// configure for Sigmoid function
	int Sigmoid_itor;
	double Sigmoid_y_inital;

	// configure for logistic regression
	double learningRate;
	int numEpochs;
	vector<double> maxItor;
	vector<double> psamples;

	//configure for HME parameters
	string p_poly_modulus;
	int p_coeff_modulus;
	int p_plain_modulus;
	vector<int> encoder_conf;

	int precision = 100;
} conf;
/*

//-------------------------< Decrease noise >----------------------------------------
Ciphertext DecreaseNoise(Ciphertext input)
{
	int buffer_length = 0;
	char *buffer = input.save(buffer_length);
  printf("!!!!!!!!!!!!!!!!1buffer length: %d\n", buffer_length);
	DecreaseNoise_SGX(eid, buffer, buffer_length);

	Ciphertext return_ans;
	return_ans.load(buffer);

	return return_ans;
}

*/

void InitialConfigure(MakeConfigure mconf)
{
	conf.Sigmoid_itor = mconf.ToInt(mconf.FindConfigure("Sigmoid_itor"));
	conf.Sigmoid_y_inital = mconf.ToDouble(mconf.FindConfigure("Sigmoid_y_inital"));
	conf.learningRate = mconf.ToDouble(mconf.FindConfigure("learningRate"));
	conf.numEpochs = mconf.ToInt(mconf.FindConfigure("runs"));
	conf.maxItor = mconf.ToDVector(mconf.FindConfigure("maxItor"));
	conf.psamples = mconf.ToDVector(mconf.FindConfigure("psamples"));
	conf.p_poly_modulus = mconf.FindConfigure("p_poly_modulus");
	conf.p_coeff_modulus = mconf.ToInt(mconf.FindConfigure("p_coeff_modulus"));
	conf.encoder_conf = mconf.ToIVector(mconf.FindConfigure("encoder_conf"));
	conf.p_plain_modulus = mconf.ToInt(mconf.FindConfigure("p_plain_modulus"));

	char* buffer = mconf.ReturnConf();
	MakeConfigure_SGX(eid, buffer, 500);

  cout<<"InitialConfigure finished."<<endl;
	//ReadData rd;
	//rd.readData(X,Y);
}

//----------------------------< Main Test Stub >---------------------------------------
int main()
{
	// Create sgx enclave
	sgx_status_t        ret = SGX_SUCCESS;
	sgx_launch_token_t  token = { 0 };
	int updated = 0;
 EncryptionParameters  parms;

	ret = sgx_create_enclave(ENCLAVE_FILE, SGX_DEBUG_FLAG, &token, &updated, &eid, NULL);
	if (ret != SGX_SUCCESS)
		return -1;

	// Initizal the Configure
	MakeConfigure mconf;
	mconf.Initalize();
	InitialConfigure(mconf);

	// Create encryption parameters
  //cout<<"***************    p_poly_modulus: "<<conf.p_poly_modulus<<"endchar"<<endl;
  parms.set_poly_modulus(conf.p_poly_modulus);
	
  parms.set_coeff_modulus(coeff_modulus_128(8192));

	parms.set_plain_modulus(conf.p_plain_modulus);
  

	// Generate keys.
	cout << "... Generating keys ..." << endl;
	generate_key_sgx(eid);
	char * public_key_buffer = new char[MAX_BUF_LEN];
	get_public_key(eid, public_key_buffer, MAX_BUF_LEN);
	public_key.load(public_key_buffer);
	cout << "... Public key generation complete ..." << endl;


#ifdef DEBUG
	//-------------< test for some function and feature >-------------------------------------------------
	cout << "This is a test for sigmoid function and sgx" << endl;
	char * secret_key_buffer = new char[MAX_BUF_LEN];
	get_secret_key(eid, secret_key_buffer, MAX_BUF_LEN);
	secret_key.load(secret_key_buffer);
	cout << "... secret key generation complete" << endl;
#endif

  SEALContext context(parms);
  
//  KeyGenerator keygen(context);
//  PublicKey public_key = keygen.public_key();
//  SecretKey secret_key = keygen.secret_key();
    
    
  Evaluator evaluator(context);
  Encryptor encryptor(context, public_key);

	IntegerEncoder encoder(context.plain_modulus());
  FractionalEncoder fracencoder(context.plain_modulus(), context.poly_modulus(),
		conf.encoder_conf[0], conf.encoder_conf[1], conf.encoder_conf[2]);
#ifdef DEBUG
	Decryptor decryptor(context, secret_key);
#endif // DEBUG

	Plaintext encLearningRate = fracencoder.encode(0.1);
  cout << "Encoded learningRate: " << encLearningRate.to_string() << endl;
  
  Ciphertext inputToSigmoid1;
  Ciphertext inputToSigmoid2;
  
  encryptor.encrypt(encoder.encode(1), inputToSigmoid1);
  encryptor.encrypt(encoder.encode(2), inputToSigmoid2);
  
  cout << "Noise budget in encrypted1: " 
        << decryptor.invariant_noise_budget(inputToSigmoid1) << " bits" << endl;
  
  for(int i = 0; i < 10; i++)
  {      
    cout<<"==================== "<<i<<" ==================== "<<endl;
    evaluator.add(inputToSigmoid1, inputToSigmoid2);
  
    cout << "Noise budget in encrypted1 + encrypted2: " 
          << decryptor.invariant_noise_budget(inputToSigmoid1) << " bits" << endl;
    evaluator.multiply(inputToSigmoid1, inputToSigmoid2);
    cout << "Noise budget in (encrypted1 + encrypted2) * encrypted2: "
          << decryptor.invariant_noise_budget(inputToSigmoid1) << " bits" << endl;
    Plaintext plain_result;
    cout << "Decrypting result: ";
    decryptor.decrypt(inputToSigmoid1, plain_result);
    cout << "Done" << endl;
    
    //Decode to obtain an integer result.
    cout << "Decoded integer: " << encoder.decode_int32(plain_result) << endl;
    
    
    // decrease noise
    int buffer_length = 0;
	  char *buffer = inputToSigmoid1.save(buffer_length);
	  DecreaseNoise_SGX(eid, buffer, buffer_length);
	  inputToSigmoid1.load(buffer);
  }

	cout << endl;
	cout << "... All is completed ..." << endl;

	// Destroy the Enclave
	if (SGX_SUCCESS != sgx_destroy_enclave(eid))
		cout << "destroy error" << endl;
	cout << "... Destroy the enclave successfully ..." << endl;

	// deallocate memory on the heap!
	delete public_key_buffer;
	delete secret_key_buffer;
  return 0;
}

