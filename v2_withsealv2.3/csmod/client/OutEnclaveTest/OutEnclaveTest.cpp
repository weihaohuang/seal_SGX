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


// socket communication
#include "../socket_client.h"

using namespace std;
using namespace seal;

#define DEBUG

PublicKey          public_key;
#if defined (DEBUG) || defined(SHOWRESULT)
SecretKey               secret_key;
#endif // DEBUG

char * secret_key_buffer;

int client_fd;

/*int connect_send() {
    while (true) {
        getline(std::cin, buf);
        if (buf == "q") {
            close(client_fd);
            break;
        } else {
            size_t size = buf.size();

            send_to_sgx(client_fd, ENCRYPT_DATA, buf.c_str(), size);

            char read_buf[1024] = {0};
            int nread = read(client_fd, read_buf, 1024);
            if (nread <= 0) {
                close(client_fd);
                break;
            } else {
                printf("recv %s.\n", read_buf);
            }

        }
    }
    return 0;
}*/


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

//	char* buffer = mconf.ReturnConf();
//	MakeConfigure_SGX(eid, buffer, 500); // ww31: be careful, I just commented out this sentence to compile client

  cout<<"InitialConfigure finished."<<endl;
	//ReadData rd;
	//rd.readData(X,Y);
}

char *int2charstring(int i)
{
  char *ret = new char[32];
  sprintf(ret, "%d", i);
  
  return ret;
}


//----------------------------< Main Test Stub >---------------------------------------
int main()
{
  EncryptionParameters  parms;
	// Initizal the Configure
	MakeConfigure mconf;
	mconf.Initalize();
	InitialConfigure(mconf);
 
  client_fd = socket_connect(IPADDRESS, PORT);  // connect to server

	parms.set_poly_modulus(conf.p_poly_modulus);
  parms.set_coeff_modulus(coeff_modulus_128(conf.p_coeff_modulus));
	parms.set_plain_modulus(conf.p_plain_modulus);
 
  SEALContext context(parms);
  
  int buffer_length = 1000;
  char *read_buf = new char[buffer_length];
  int nread;

   
  // send encryption parameters
  send_to_sgx(client_fd, ENC_PARAMETER_POLYMOD, parms.poly_modulus().to_string().c_str(), parms.poly_modulus().to_string().length());
  nread = read(client_fd, read_buf, buffer_length);
  if (nread <= 0) {
      close(client_fd);
      return -1;
  } else {
      printf("recv %s.\n", read_buf);
      printf("<<<<<<<<<<<<<<<  recv chars %d\n", read_buf[0]);
  }

  send_to_sgx(client_fd, ENC_PARAMETER_COEFMOD, int2charstring(8192), strlen(int2charstring(8192)));
  nread = read(client_fd, read_buf, buffer_length);
  if (nread <= 0) {
      close(client_fd);
      return -1;
  } else {
      printf("recv %s.\n", read_buf);
      printf("<<<<<<<<<<<<<<<  recv chars %d %d\n", read_buf[0], read_buf[100]);
  }
  send_to_sgx(client_fd, ENC_PARAMETER_PLAINMOD, int2charstring(conf.p_plain_modulus), strlen(int2charstring(conf.p_plain_modulus)));
  nread = read(client_fd, read_buf, buffer_length);
  if (nread <= 0) {
      close(client_fd);
      return -1;
  } else {
      printf("recv %s.\n", read_buf);
      printf("<<<<<<<<<<<<<<<  recv chars %d %d\n", read_buf[0], read_buf[100]);
  }
   
//  cout << "parms.poly_modulus() count: " << parms.poly_modulus().coeff_count() << endl;
//  cout << "parms.poly_modulus(): " << parms.poly_modulus().to_string() << endl;
//	cout << "parms.coeff_modulus(): " << parms.coeff_modulus().to_string() << endl;
//	cout << "parms.plain_modulus(): " << parms.plain_modulus().to_string() << endl;
	// Generate keys.
	cout << "... Generating keys ..." << endl;
	KeyGenerator generator(context);
 
  public_key = generator.public_key();
	cout << "... Public key generation complete ..." << endl;

  // send public key
  int buffer_size = 0;
	char* tmp_p_k_b = public_key.save(buffer_size);
  send_to_sgx(client_fd, PUBLIC_KEY, (const char *)tmp_p_k_b, buffer_size);
  nread = read(client_fd, read_buf, buffer_length);
  if (nread <= 0) {
      close(client_fd);
      return -1;
  } else {
      printf("recv %s.\n", read_buf);
      printf("<<<<<<<<<<<<<<<  recv chars %d %d\n", read_buf[0], read_buf[100]);
  }
  
#ifdef DEBUG
	//-------------< test for some function and feature >-------------------------------------------------
  secret_key = generator.secret_key();
	
	cout << "... secret key generation complete" << endl;
  // send private key
  buffer_size = 0;
	char* tmp_s_k_b = secret_key.save(buffer_size);
  send_to_sgx(client_fd, PRIVATE_KEY, (const char *)tmp_s_k_b, buffer_size);
  nread = read(client_fd, read_buf, buffer_length);
  if (nread <= 0) {
      close(client_fd);
      return -1;
  } else {
      printf("recv %s.\n", read_buf);
      printf("<<<<<<<<<<<<<<<  recv chars %d %d\n", read_buf[0], read_buf[100]);
  }
#endif
  Evaluator evaluator(context);
  Encryptor encryptor(context, public_key);
  
  IntegerEncoder encoder(context.plain_modulus());
  FractionalEncoder fracencoder(context.plain_modulus(), context.poly_modulus(),
		conf.encoder_conf[0], conf.encoder_conf[1], conf.encoder_conf[2]);
   
#ifdef DEBUG
	Decryptor decryptor(context, secret_key);
#endif // DEBUG

	Plaintext encLearningRate = fracencoder.encode(0.1);
  cout << "A simple test, encoded learningRate: " << encLearningRate.to_string() << endl;

  Ciphertext inputToSigmoid1;
  Ciphertext inputToSigmoid2;
  
  encryptor.encrypt(encoder.encode(1), inputToSigmoid1);
  encryptor.encrypt(encoder.encode(2), inputToSigmoid2);
 
  SimulationEvaluator sim_evaluator;
  
  Simulation sim_input1 = sim_evaluator.get_fresh(parms, conf.p_coeff_modulus, conf.p_coeff_modulus);
  Simulation sim_input2 = sim_evaluator.get_fresh(parms, conf.p_coeff_modulus, conf.p_coeff_modulus);
  
  
  cout << "Noise in the inputToSigmoid: " << decryptor.invariant_noise_budget(inputToSigmoid1)<<" bits" << endl;
  cout<<"Estimted noise (begin): "<<sim_input1.invariant_noise_budget()<<endl;
	int i;
 
	for (i = 0; i < 10; i++)
	{     
		evaluator.add(inputToSigmoid1, inputToSigmoid2);
		evaluator.multiply(inputToSigmoid1, inputToSigmoid2);
   
    sim_input1 = sim_evaluator.add(sim_input1, sim_input2);
    sim_input1 = sim_evaluator.multiply(sim_input1, sim_input2);
    
    cout << "Estimated noise: " << sim_input1.invariant_noise_budget() << " bits" << endl;
    
    char message[32] = {0};
    sprintf(message, "%d", sim_input1.invariant_noise_budget());
    send_to_sgx(client_fd, STATUS_NOISE, message, 32);
    
    char response[32] = {0};
    
    int mread = read(client_fd, response, 32);
    if (mread <= 0) {
      close(client_fd);
      return -1;
    } else {
      printf("recv %s.\n", response);
    }
     
    int buffer_length = 0;
	  char *buffer = inputToSigmoid1.save(buffer_length);
    
//  	printf("<<<<<<<<<<<<<<<  send chars %d %d\n", buffer[0], buffer[100]);
    send_to_sgx(client_fd, ENCRYPT_DATA, buffer, buffer_length);
    
    sim_input1 = sim_evaluator.get_fresh(parms, conf.p_coeff_modulus, conf.p_coeff_modulus);
    cout << "Estimated noise: " << sim_input1.invariant_noise_budget()<<" bits" << endl;
    
    char *read_buf = new char[buffer_length];
    char *tmp = read_buf;
    int total_length = 0;
recv:
    int nread = read(client_fd, tmp, buffer_length);
    
    if (nread <= 0) {
        close(client_fd);
        break;
    } else {
        total_length += nread;
        if(total_length < buffer_length) 
        {
          tmp += nread;
          goto recv;
        }
//        printf("recv %s.\n", read_buf);
//        printf("<<<<<<<<<<<<<<<  recv chars %d %d, length: %d\n", read_buf[0], read_buf[100], total_length);
        inputToSigmoid1.load(read_buf);
    }
    delete [] read_buf;
   
#ifdef DEBUG

		Plaintext ans;
    decryptor.decrypt(inputToSigmoid1, ans);

		
		cout << "Noise in the inputToSigmoid: " << decryptor.invariant_noise_budget(inputToSigmoid1)<<" bits" << endl;
#endif
	}

	cout << endl;
	cout << "... All is completed ..." << endl;

	// deallocate memory on the heap!
#ifdef DEBUG
  delete [] secret_key_buffer;
  delete [] tmp_s_k_b;
#endif
  delete [] tmp_p_k_b;
  delete [] read_buf;
  close(client_fd);
  return 0;
}




