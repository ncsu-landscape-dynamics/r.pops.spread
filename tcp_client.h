#include<string>  //string
#include<sys/socket.h>    //socket
#include<arpa/inet.h> //inet_addr

using std::string;

/**
    TCP Client class
*/
class tcp_client
{
private:
    int sock;
    string address;
    int port;
    struct sockaddr_in server;
    
public:
    tcp_client();
    bool conn(string, int);
    bool send_data(string data);
    string receive(int, int &error);
};

