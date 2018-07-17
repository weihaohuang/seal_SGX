#include "SealEnclaveTest_t.h"

#include "sgx_trts.h"

#include <sstream>
//#include <cstring>
#include <chrono>
#include <vector>
#include "../seal/seal.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cmath>
#include <map>

using namespace std;
using namespace seal;

//------------------------------------------------------------------------------
/*
* printf:
*   Invokes OCALL to display the enclave buffer to the terminal.
*/
//------------------------------------------------------------------------------
void printf(const char *fmt, ...) {
	char buf[BUFSIZ] = { '\0' };
	va_list ap;
	va_start(ap, fmt);
	vsnprintf(buf, BUFSIZ, fmt, ap);
	va_end(ap);
	ocall_print(buf);
}

// SEALContext context(parms);

struct Configure
{
	string p_poly_modulus;
	int p_coeff_modulus;
	int p_plain_modulus;
};


map<int, struct Configure> conf;
//EncryptionParameters parms_sgx;
map<int, SecretKey> secret_key_sgx;
//BigPoly secret_key_sgx;
map<int, PublicKey> public_key;
//BigPolyArray public_key;



double decrypted_number;

int check_Index()
{
	int flag = 0;
	if (decrypted_number > 0)
		flag = 0;
	else
		flag = 1;
	return flag;
}

void sigmod_sgx(int client_id, char* buffer, size_t len,int trainingSize,int precision)
{	
}

void DecreaseNoise_SGX(int client_id, char* buf, size_t len)
{
  // test whether the public/private keys have been set
  // if yes, retrieve the keys; otherwise, print error message
  EncryptionParameters parms_sgx;

  parms_sgx.set_poly_modulus(conf[client_id].p_poly_modulus);
  parms_sgx.set_coeff_modulus(coeff_modulus_128(conf[client_id].p_coeff_modulus));
	parms_sgx.set_plain_modulus(conf[client_id].p_plain_modulus);
  
  SEALContext context(parms_sgx);
  
	Encryptor encryptor(context, public_key[client_id]);
	Decryptor decryptor(context, secret_key_sgx[client_id]);

//	Evaluator evaluator(parms_sgx);
	Plaintext encoded_number;
	Ciphertext encrypted_rational;
	Plaintext plain_result;

	encrypted_rational.load(buf);
	decryptor.decrypt(encrypted_rational, plain_result);
 
  // ww31: it may need to decode & encode, however I remove it for now
	encryptor.encrypt(plain_result, encrypted_rational);

	int length = 0;
	char* tmp_buf = encrypted_rational.save(length);

	memcpy(buf, tmp_buf, length);

	delete[] tmp_buf;

}

// ********************** CHECK WITH YONGSOO ON EVALUATING THE ~0,1 VALS WITHIN plainModBound space *****************
void AddInRow_SGX(char* buf, size_t len,int trainingSize,int precision)
{
}

void set_public_key(int client_id, char* public_key_buffer, size_t len)
{
  public_key[client_id].load(public_key_buffer);
}

void set_secret_key(int client_id, char* secret_key_buffer, size_t len)
{
  secret_key_sgx[client_id].load(secret_key_buffer);
}


void MakeConfigure_SGX(int client_id, char* polymod, int polymodlen, char* coefmod, int coefmodlen, char* plainmod, int plainmodlen)
{
  EncryptionParameters  parms_sgx;
  
  polymod[polymodlen] = 0;
  coefmod[coefmodlen] = 0;
  plainmod[plainmodlen] = 0;
  printf("................. client id: %d, polymod: %s, length: %d\n", client_id, polymod, polymodlen);
  
//  parms_sgx.poly_modulus() = polymod;
 // parms_sgx.coeff_modulus() = coefmod;
 // parms_sgx.plain_modulus() = atoi(plainmod);
 
 	conf[client_id].p_poly_modulus = polymod;
  conf[client_id].p_coeff_modulus = atoi(coefmod);
  conf[client_id].p_plain_modulus = atoi(plainmod);
  
//	memcpy(&parms_sgx, ConfigureBuffer, len);
//  printf("sizeof(EncryptionParameters): %d, len: %d\n", sizeof(EncryptionParameters), len);
  printf("parms_sgx.coeff_modulus_: %s\n", coefmod);
  printf("parms_sgx.plain_modulus_: %d\n", atoi(plainmod));
}
