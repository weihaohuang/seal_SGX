This project implements a client (who performs homomorphic add & multiplication) and a server (who performs bootstrapping).  

Whenever a bootstrapping is needed (i.e. the noise exceeds a threshold, which can be estimated using the evaluation key), the ciphertext is sent to the server. The server will decrypt it and re-encrypt it (equivalent to a bootstrapping), and send it to the client. The project uses socket communication between the client and server.   

As the first step before the communication, the server needs to get the homomorphic encryption key configuration, and the private key (after a remote attestation is performed and the key agreement is established).

### Status
- [X] Communication established.
- [X] Data transfer success.
- [X] Encryption configuration transfer success.
- [X] Decrease_noise   
After ``decrease_noise`` when the ciphertext is transferred back to the client, the client cannot load it successfully.(_fixed_) 
- [X] Ciphertext after noise decrease sent back to client
- [X] Multiple clients support
  - Needs to bind the keys to each client
  
### Missing
- [ ] Multiple server enclaves
  - which may have a negtive effect when multiple threads run in a same enclave
- [ ] Scheduling multiple client requests
  - A simple scheduling method
    - The server maintains a task queue with priorities
    - The client sents the current distance to the threshold after each (or several) homomorphic computation(s)
    - The server decides which client can send the ciphertext for bootstrapping: send the decision to the client
  - Take relinearization (decreasing the size of the ciphertext and slightly increasing the noise), or Galois Automorphisms (supporting batching operations) into consideration to minimize network overhead
- [ ] Removing side channel leakages
  - Understand current leakage
  - No secret dependent branches
  - No secret dependent memory accesses (e.g. secret as array index)
  - Memory pool
- [ ] Evaluation
  - Types
    - bootstrapping
    - single client-server task v.s. leveled HME with bigger parameters sizes
    - multiple client-server tasks: evaluating scheduling algorithm
  - Evaluation tasks
    - Logistic regression
- [ ] Bug fixes
  - Even if ``decrease_noise`` is not called, the client & server channel sometimes terminate unexpectly.
  - Possible crashes when the data trasfered through socket communication is incorrect. (rare cases)   
    Temporarily fixed by ingnoring buffer length is < 0 or > length

### Future work
- Support more HME schemes, at least one GPU-based HME scheme implementation
