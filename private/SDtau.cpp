#include "VineCPP_header.hpp"
#include <matrix.h>
#include <mex.h>

// Helper function: Compare two pairs of doubles
bool UV_sorter(std::pair<double, double> const& lhs, std::pair<double, double> const& rhs) {
    if (lhs.first != rhs.first)
        return lhs.first < rhs.first;
    return lhs.second < rhs.second;
}

// Helper function: Sort two vectors; primarily with respect to the first and secondarily with respect to the second
void SortUV(std::pair<double, double> UV[], double *U, double *V,unsigned int n)
{
    unsigned int i;
    
    for (i=0;i<n;i++)
    {
        UV[i] = std::make_pair(U[i],V[i]);
    }
    
    std::sort(&UV[0],&UV[n], &UV_sorter);
    return;
}

// AVL tree node struct and its constructor
struct AVL_Node{
    double Value; // The value of the node
    struct AVL_Node *Left; // Node to the left
    struct AVL_Node *Right; // Node to the right
    unsigned int Height; // The height of the tree at the node
    unsigned int Size; // The size of the tree at the node
    unsigned int Count; // The number of times the node appears
    AVL_Node(double V)
    {
        Value = V;
        Left = NULL;
        Right = NULL;
        Height = 1;
        Size = 1;
        Count = 1;
    }
};

void Print_AVL_Node(AVL_Node *Node){
    printf("Node Value: %2.2f; Height: %i; Size: %i; Count: %i \n",Node->Value,Node->Height,Node->Size,Node->Count);
}

// Helper functions to obtain the height or size
unsigned int height(AVL_Node *Node)
{
    if (Node==0)
        return 0;
    return Node->Height;
}

unsigned int size(AVL_Node *Node)
{
    if (Node==0)
        return 0;
    return Node->Size;
}

// Helper function to compute the balance
int balance(AVL_Node *Node)
{
    if (Node==0)
        return 0;
    return height(Node->Left) - height(Node->Right);
}


// Helper function to rotate the subtree at node Node1 to the right
AVL_Node *Rotate_Right(AVL_Node *Node1)
{
    // Getting the tree to the left and its right one
    AVL_Node *Node2 = Node1->Left;
    AVL_Node *Node3 = Node2->Right;
    // Rotation
    Node2->Right = Node1;
    Node1->Left = Node3;
    
    // Compute the new height and size values
    Node1->Height = 1 + max(height(Node1->Left),height(Node1->Right));
    Node1->Size = Node1->Count + size(Node1->Left) + size(Node1->Right);
            
    Node2->Height = 1 + max(height(Node2->Left),height(Node2->Right));
    Node2->Size = Node2->Count + size(Node2->Left) + size(Node2->Right);
    
    return Node2;
}

// Helper function to rotate the subtree at node Node1 to the left
AVL_Node *Rotate_Left(AVL_Node *Node1)
{
    // Getting the tree to the right and its left one
    AVL_Node *Node2 = Node1->Right;
    AVL_Node *Node3 = Node2->Left;
    // Rotation
    Node2->Left = Node1;
    Node1->Right = Node3;
    
    // Compute the new height and size values
    Node1->Height = 1 + max(height(Node1->Left),height(Node1->Right));
    Node1->Size = Node1->Count + size(Node1->Left) + size(Node1->Right);
            
    Node2->Height = 1 + max(height(Node2->Left),height(Node2->Right));
    Node2->Size = Node2->Count + size(Node2->Left) + size(Node2->Right);
    
    return Node2;
}

// Insert a new node into the AVL tree
AVL_Node* AVL_insert(AVL_Node* Node, double V, unsigned int *NumbBefore, unsigned int *NumbEqual)
{
    int Balance;
    
    if (Node==NULL)
    {
        AVL_Node *N = new AVL_Node(V);
        return(N);
    }
    
    if (V == Node->Value)
    {
        Node->Count = Node->Count +1;
            
        *NumbEqual =  Node->Count;
        *NumbBefore = *NumbBefore + size(Node->Left);
    }
    else if (V < Node->Value)
    {
        Node->Left = AVL_insert(Node->Left,V,NumbBefore,NumbEqual);
    }
    else
    {
        Node->Right = AVL_insert(Node->Right,V,NumbBefore,NumbEqual);
        
        *NumbBefore = *NumbBefore + size(Node->Left) + Node->Count;
    }
    
    Node->Height = 1 + max(height(Node->Left),height(Node->Right));
    Node->Size = Node->Count + size(Node->Left) + size(Node->Right);
    
    Balance = balance(Node);
    
    if (Balance>1 && V < Node->Left->Value)
        return Rotate_Right(Node);
    
    if (Balance<-1 && V > Node->Right->Value)
        return Rotate_Left(Node);
    
    if (Balance>1 && V > Node->Left->Value)
    {
        Node->Left = Rotate_Left(Node->Left);
        return Rotate_Right(Node);
    }
    
    if (Balance<-1 && V < Node->Right->Value)
    {
        Node->Right = Rotate_Right(Node->Right);
        return Rotate_Left(Node);
    }
    
    return Node;
}

void AVL_Tree (std::pair<double, double> UV[], unsigned int *NumbBefore, unsigned int *NumbEqual, unsigned int n)
{
    unsigned int i;
    AVL_Node *AVL_Root_Node = NULL;
    
    for (i=0;i<n;i++)
    {
        NumbBefore[i] = 0;
        NumbEqual[i] = 1;
        AVL_Root_Node = AVL_insert(AVL_Root_Node,UV[i].second,&NumbBefore[i],&NumbEqual[i]);
        
    }
    
    delete [] AVL_Root_Node;
    
}
void SD_Kendall_Tau_Matrix(double *tau, double *U, unsigned int d, unsigned int n)
{
    unsigned int i,j;
    
    for (i=0;i<d;i++)
    {
        tau[i*d+i] = 1;
        for (j=0;j<i;j++)
        {
            tau[i*d+j] = SD_Kendall_Tau(&U[i*n],&U[j*n],n);
            tau[j*d+i] = tau[i*d+j];
        }
    }
    return;
}

double SD_Kendall_Tau(double *U, double *V, unsigned int n)
{
    double tau=0;
    
    unsigned int i, *NumbBefore, *NumbEqual;
    unsigned int NBefore[n-1], NEqual[n-1];
    NumbBefore = &NBefore[0];
    NumbEqual = &NEqual[0];
    std::pair<double, double> UV[ n ];
    
    SortUV(UV,U,V,n);
    AVL_Tree(UV,NumbBefore,NumbEqual,n);
    
    int Concordant = 0, Discordant = 0, ExtraX = 0, ExtraY = 0, ACount = 0, BCount = 0, CCount = 0, DCount = 0, ECount = 0;
    double h1=0, h2=0;
    
    double PrevU = UV[0].first-1;
    double PrevV = UV[0].second-1;
    
    
    for (i=0;i<n;i++)
    {
        if (UV[i].first != PrevU)
        {
            DCount = 0;
            ECount = 1;
        }
        else
        {
            if (UV[i].second == PrevV)
            {
                ECount++;
            }
            else
            {
                DCount += ECount;
                ECount = 1;
            }
        }
        
        ACount = NumbBefore[i] - DCount;
        BCount = NumbEqual[i] - ECount;
        CCount = i - (ACount + BCount + DCount + ECount - 1);
        
        ExtraY += DCount;
        ExtraX += BCount;
        
        Concordant += ACount;
        Discordant += CCount;
        
        PrevU = UV[i].first;
        PrevV = UV[i].second;
    }
    
    h1 = Concordant + Discordant + ExtraX;
    h2 = Concordant + Discordant + ExtraY;
    double h3 = Concordant - Discordant;
    
    tau = h3/sqrt(h1*h2);
    
    return tau;
}


