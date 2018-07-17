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


SecretKey secret_key_sgx;
PublicKey public_key;
double decrypted_number;

static struct configure_SGX
{
	string p_poly_modulus;
	int p_coeff_modulus;
	vector<int> encoder_conf;
	int p_plain_modulus;
}conf_SGX;

void DecreaseNoise_SGX(char* buf, size_t len)
{
  EncryptionParameters parms_sgx;
  parms_sgx.set_poly_modulus(conf_SGX.p_poly_modulus);
  parms_sgx.set_coeff_modulus(coeff_modulus_128(8192));
	parms_sgx.set_plain_modulus(conf_SGX.p_plain_modulus);
 
  SEALContext context(parms_sgx);
  
	Encryptor encryptor(context, public_key);
	Decryptor decryptor(context, secret_key_sgx);

	Ciphertext encrypted_rational;
	Plaintext plain_result;

	encrypted_rational.load(buf);
	decryptor.decrypt(encrypted_rational, plain_result);
	
	encryptor.encrypt(plain_result, encrypted_rational);

	int length = 0;
	char* tmp_buf = encrypted_rational.save(length);

	memcpy(buf, tmp_buf, length);

	delete[] tmp_buf;

}


//----------------------------< operation for keys >-----------------------------------------
void generate_key_sgx()
{
  EncryptionParameters parms_sgx;
  parms_sgx.set_poly_modulus(conf_SGX.p_poly_modulus);
  parms_sgx.set_coeff_modulus(coeff_modulus_128(8192));
	parms_sgx.set_plain_modulus(conf_SGX.p_plain_modulus);

  
  SEALContext context(parms_sgx);

	//generate the public_key and secret_key;
	KeyGenerator generator(context);
	//generator.generate();
 
  PublicKey tmp_public_key = generator.public_key();
  SecretKey tmp_secret_key_sgx = generator.secret_key();
    
	//BigPolyArray tmp_public_key = generator.public_key();
	//BigPoly tmp_secret_key_sgx = generator.secret_key();

	// Store public_key and secret_key in enclave
	int buffer_size = 0;
	char* tmp_s_k_b = tmp_secret_key_sgx.save(buffer_size);
	secret_key_sgx.load(tmp_s_k_b);

	buffer_size = 0;
	char* tmp_p_k_b = tmp_public_key.save(buffer_size);
	public_key.load(tmp_p_k_b);

	delete[] tmp_s_k_b;
	delete[] tmp_p_k_b;
}

void get_public_key(char* public_key_buffer,size_t len)
{
	int buffer_size = 0;
	char* tmp_p_k_b = public_key.save(buffer_size);
	memcpy(public_key_buffer, tmp_p_k_b, buffer_size);
  
	delete[] tmp_p_k_b;
}

void get_secret_key(char* secret_key_buffer,size_t len)
{
	int buffer_size = 0;
	char* tmp_s_k_b = secret_key_sgx.save(buffer_size);
	memcpy(secret_key_buffer, tmp_s_k_b, buffer_size);
	delete[] tmp_s_k_b;
}


// Initial the Configure in the SGX
string FindConfigure(string input,char* ConfigureBuffer)
{
	string ans;
	string tCon = ConfigureBuffer;
	string temp;

	for (int i = 0; i < 500; i++)
	{
		if (ConfigureBuffer[i] == input[0])
		{
			temp = tCon.substr(i, input.length());
			if (temp == input)
			{
				int count = 1;
				i += input.length();
				while (ConfigureBuffer[i] != '=')
					i++;
				while (ConfigureBuffer[i + count] != '#')
					count++;
				ans = tCon.substr(i + 1, count - 1);
				break;
			}
		}
		else
		{
			while (ConfigureBuffer[i] != '#')
				i++;
		}
	}

	return ans;
}

int ToInt(string input)
{
	int value = atoi(input.c_str());
	return value;
}

vector<int> ToIVector(string input)
{
	vector<int> ans;
	int position = 0;
	for (int i = 0; i < input.length(); i++)
	{
		if (input[i] == ';')
		{
			string temp = input.substr(position, i - position);
			ans.push_back(ToInt(temp));
			position = i + 1;
		}
		if (i == input.length() - 1)
		{
			string temp = input.substr(position, i - position + 1);
			ans.push_back(ToInt(temp));
		}
	}
	return ans;
}

void MakeConfigure_SGX(char* ConfigureBuffer, size_t len)
{
	conf_SGX.encoder_conf=ToIVector(FindConfigure("encoder_conf",ConfigureBuffer));
	conf_SGX.p_poly_modulus = FindConfigure("p_poly_modulus", ConfigureBuffer);
	conf_SGX.p_coeff_modulus = ToInt(FindConfigure("p_coeff_modulus", ConfigureBuffer));
	conf_SGX.p_plain_modulus = ToInt(FindConfigure("p_plain_modulus", ConfigureBuffer));
}
