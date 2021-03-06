#include "SealEnclaveTest_u.h"
#include "../socket_server.h"
#include "sgx_urts.h"
#include <stdio.h>
#include <string.h>
#include "Matrix.h"
#include <vector>
#include <algorithm>
#include <math.h>
#include <iostream>
#include <sstream>
#include <chrono>
#include <unordered_map>
#include <assert.h>
#include "TestData.h"
#include <map>

using namespace std;

#define ENCLAVE_FILE "SealEnclaveTest.signed.so"
#define MAX_BUF_LEN 600000
#define round(x) (x - floor(x) >= 0.5 ? floor(x) + 1 : floor(x)) 
#define SIGNIFICANT_FIGURES 2
#define Simplify_NewSigmoid
#define DEBUG
#define SHOWRESULT

//------------------------------------------------------------------------------
void ocall_print(const char* str) {
	printf("\033[96m%s\033[0m", str);
}

//vector<vector<double>> X;
//vector<double> Y;

//EncryptionParameters  parms;
sgx_enclave_id_t      eid;

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

//---------< Round Function >----------------------------------------------------------
double Round(double value)
{
	int n = SIGNIFICANT_FIGURES;
	return round(value * pow(10, n)) / pow(10, n);
}

//---------< Create vector by step >--------------------------------------------------
void Create_Vector_By_Step(vector<double>& tmpVector, double start, double step, double end)
{
	for (double i = start; i <= end; i += step)
		tmpVector.push_back(i);
}

//-------------------------< Decrease noise >----------------------------------------
/*BigPolyArray DecreaseNoise(BigPolyArray input)
{
	int buffer_length = 0;
	char *buffer = input.save(buffer_length);
  printf("!!!!!!!!!!!!!!!!1buffer length: %d\n", buffer_length);
	DecreaseNoise_SGX(eid, buffer, buffer_length);

	BigPolyArray return_ans;
	return_ans.load(buffer);

	return return_ans;
}*/

//----------< Random Function >----------------------------------------------------------
double GetRandom(double min, double max) {
	/* Returns a random double between min and max */
	return ((double)rand()*(max - min) / (double)RAND_MAX - min);
}


//-------------< Addition in SGX >------------------------------------------------------
/*BigPolyArray AddInRow(BigPolyArray input,int trainingSize)
{
	int buffer_length = 0;
	char *buffer = input.save(buffer_length);

	AddInRow_SGX(eid, buffer, buffer_length, trainingSize,conf.precision);

	BigPolyArray return_ans;
	return_ans.load(buffer);

	delete[] buffer;

	return return_ans;
}*/

//-----------< sigmod function based on hme >-------------------------------------------
/*BigPolyArray sigmod_Hme(BigPolyArray input,int trainingSize)
{
	int buffer_length = 0;
	char *buffer = input.save(buffer_length);

	sigmod_sgx(eid, buffer, buffer_length, trainingSize,conf.precision);

	BigPolyArray return_ans;
	return_ans.load(buffer);

	return return_ans;
}*/

int main()
{
  // Create sgx enclave
  sgx_status_t        ret = SGX_SUCCESS;
  sgx_launch_token_t  token = { 0 };
  int updated = 0;
  ret = sgx_create_enclave(ENCLAVE_FILE, SGX_DEBUG_FLAG, &token, &updated, &eid, NULL);
  if (ret != SGX_SUCCESS)
  	return -1;
  cout<<"Enclave loaded."<<endl;
  // Initizal the Configure
  
  //	parms.poly_modulus() = conf.p_poly_modulus;
  //	parms.coeff_modulus() = ChooserEvaluator::default_parameter_options().at(conf.p_coeff_modulus);
  //	parms.plain_modulus() = conf.p_plain_modulus;
  int listen_fd = socket_bind(IPADDRESS, PORT);
  int max_fd = -1;
  int nready;
  fd_set readfds;
  int clients_fd[IPC_MAX_CONN];
  
  memset(clients_fd, -1, sizeof(clients_fd));
  
  while (true) {
  
      FD_ZERO(&readfds);
      FD_SET(listen_fd, &readfds);
      max_fd = listen_fd;
  
      for (size_t i=0; i < IPC_MAX_CONN; i++)
      {
          if (clients_fd[i] != -1) {
              FD_SET(clients_fd[i], &readfds);
              max_fd = clients_fd[i] > max_fd ? clients_fd[i] : max_fd;
          }
      }
      nready = select(max_fd+1, &readfds, NULL, NULL, NULL);
      if (nready == -1) {
          perror("select error.");
          return 1;
      }
      if (FD_ISSET(listen_fd, &readfds)) {
          accpet_client(clients_fd, listen_fd);
      }
      recv_client_msg(clients_fd, &readfds);
  }
  
  // Destroy the Enclave
  if (SGX_SUCCESS != sgx_destroy_enclave(eid))
  	cout << "destroy error" << endl;
  cout << "... Destroy the enclave successfully ..." << endl;
  
  
  return 0;
}

// multiple clients
bool establish_connection[IPC_MAX_CONN] = {false};
map<int, bool> flag_polymod, flag_coefmod, flag_plainmod, configured;
map<int, int> polymodlen, coefmodlen, plainmodlen;
map<int, char *> polymod, coefmod, plainmod; 


/* single client  
//bool flag_polymod = false, flag_coefmod = false, flag_plainmod = false;
//char polymod[100] = {0}, coefmod[100] = {0}, plainmod[100] = {0}; 
//int polymodlen = 0, coefmodlen = 0, plainmodlen = 0;
//bool configured = false;
*/

// message processing
/**
 * 通过select查询到fdset之后,循环遍历每个fd是否就绪
 * @param clients_fd
 * @param readfds
 */
void recv_client_msg(int *clients_fd, fd_set *readfds) {
    char *buf = new char[bufflen];
    struct message_head head;

    for (size_t i = 0; i < IPC_MAX_CONN; ++i) {
        if (clients_fd[i] == -1) {
            continue;
        } else if (FD_ISSET(clients_fd[i], readfds)) {
            if(!establish_connection[i])
            {
              establish_connection[i] = true;
              printf("client id: %d\t%d\n", i, clients_fd[i]);
              flag_polymod[i] = flag_coefmod[i] = flag_plainmod[i] = false; // set the flag to false when a new client is connected.
              polymodlen[i] = 0, coefmodlen[i] = 0, plainmodlen[i] = 0;
              polymod[i] = new char[100];
              coefmod[i] = new char[100];
              plainmod[i] = new char[100];
            }
            int n = read(clients_fd[i], &head, sizeof(struct message_head));
            if (n <= 0) {
                FD_CLR(clients_fd[i], readfds);
                printf("one socket close\n");
                close(clients_fd[i]);
                clients_fd[i] = -1;
                establish_connection[i] = false;
                
                delete [] polymod[i];
                delete [] coefmod[i];
                delete [] plainmod[i];
              
                continue;
            }
            printf("command is: %d, buffer length: %d\n", head.cmd, head.data_len);
            if(head.data_len > bufflen || head.data_len < 0)
            {
              printf("overflow detected.\n");
              continue;
            }
            read(clients_fd[i], buf, head.data_len);
            if(head.cmd == ENC_PARAMETER_POLYMOD)
            {
              flag_polymod[i] = true;
              strncpy(polymod[i], buf, head.data_len);
              polymodlen[i] = head.data_len;
              printf(">>>>>>>>>>>>>>   recv polymod: %s, %s, length: %d\n", polymod[i], buf, polymodlen[i]);
              if(!configured[i] && flag_polymod[i] && flag_coefmod[i] && flag_plainmod[i])
              {
                configured[i] = true;
                MakeConfigure_SGX(eid, clients_fd[i], polymod[i], polymodlen[i], coefmod[i], coefmodlen[i], plainmod[i], plainmodlen[i]);
              }
            }else if(head.cmd == ENC_PARAMETER_COEFMOD)
            {
              flag_coefmod[i] = true;
              strncpy(coefmod[i], buf, head.data_len);
              coefmodlen[i] = head.data_len;
              if(!configured[i] && flag_polymod[i] && flag_coefmod[i] && flag_plainmod[i])
              {
                configured[i] = true;
                MakeConfigure_SGX(eid, clients_fd[i], polymod[i], polymodlen[i], coefmod[i], coefmodlen[i], plainmod[i], plainmodlen[i]);
              }
            }else if(head.cmd == ENC_PARAMETER_PLAINMOD)
            {
              flag_plainmod[i] = true;

              strncpy(plainmod[i], buf, head.data_len);
              plainmodlen[i] = head.data_len;
              if(!configured[i] && flag_polymod[i] && flag_coefmod[i] && flag_plainmod[i])
              {
                configured[i] = true;
                MakeConfigure_SGX(eid, clients_fd[i], polymod[i], polymodlen[i], coefmod[i], coefmodlen[i], plainmod[i], plainmodlen[i]);
              }
            }
            else if(head.cmd == PUBLIC_KEY)
            {
              set_public_key(eid, clients_fd[i], buf, head.data_len);
            }else if(head.cmd == PRIVATE_KEY)
            {
              set_secret_key(eid, clients_fd[i], buf, head.data_len);
            }else if(head.cmd == ENCRYPT_DATA)
            {
              printf("<<<<<<<<<<<<  Decrease noise in SGX enclave called....<<<<<<last char: %d %d\n", buf[0], buf[100]);
              DecreaseNoise_SGX(eid, clients_fd[i], buf, head.data_len);
            }else if(head.cmd == STATUS_NOISE)
            {
              printf("current noise budget: %d\n", atoi(buf));
              // decide whether to request for bootstrapping
            }else
            {
              printf(">>> unknown command.\n");
            }
            
            if(head.cmd == ENCRYPT_DATA)
            {
              write(clients_fd[i], buf, head.data_len);
              printf("<<<<<<<<<<<<  Written back ciphertext to client.<<<<<<last char: %d %d\n", buf[0], buf[100]);
            }else
            {
              char *ret = "processed";
              write(clients_fd[i], ret, strlen(ret));
              printf("<<<<<<<<<<<<  Written back processed to client.");
              //handle_client_msg(clients_fd[i], &head, buf);
            }
            printf("loop finished.\n");
        }
    }
    delete [] buf;
}