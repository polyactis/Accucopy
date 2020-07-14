/*
 * 2013.09.21 Yu Huang, copyright. a red-black tree template, c++
 * based on http://web.mit.edu/~emin/www.old/source_code/cpp_trees/index.html (check c++/trees folder in test repository)
 *
 * asoliman's version https://code.google.com/p/rbtrees/ is too buggy.
 * so try luck on this version.
 *
 * Examples
 *
typedef set<BreakPoint*> rbNodeDataType;
typedef RedBlackTree<BreakPointKey, rbNodeDataType > treeType;
typedef RedBlackTreeNode<BreakPointKey, rbNodeDataType > rbNodeType;

	treeType rbTree = treeType();
	BreakPoint* leftBreakPointPtr=NULL;
	BreakPoint* rightBreakPointPtr=NULL;

	rbNodeType* currentNodePtr=rbTree.nil;
	...
	currentNodePtr = rbTree.queryTree(bpKey);
	if (rbTree.isNULLNode(currentNodePtr)){
		rbNodeDataType* dataPtr = new rbNodeDataType();
		currentNodePtr = rbTree.insertNode(bpKey, dataPtr);
	}
	currentNodePtr->getDataPtr()->insert(bpPtr);
 */

#ifndef E_REDBLACK_TREE
#define E_REDBLACK_TREE

#include<math.h>
#include<limits.h>
#include<iostream>
#include <boost/format.hpp>
#include "misc.h"
#include <cstring>      //for strcat
using namespace std;

//  CONVENTIONS:
//                Function names: Each word in a function name begins with
//                a capital letter.  An example function name is
//                CreateRedTree(a,b,c). Furthermore, each function name
//                should begin with a capital letter to easily distinguish
//                them from variables.
//
//                Variable names: Each word in a variable name begins with
//                a capital letter EXCEPT the first letter of the variable
//                name.  For example, int newLongInt.  Global variables have
//                names beginning with "g".  An example of a global
//                variable name is gNewtonsConstant.

#ifndef MAX_INT
#define MAX_INT INT_MAX // some architechturs define INT_MAX not MAX_INT
#endif
const int MIN_INT = -MAX_INT;

// The RedBlackEntry class is an Abstract Base Class.  This means that no
// instance of the RedBlackEntry class can exist.  Only classes which
// inherit from the RedBlackEntry class can exist.  Furthermore any class
// which inherits from the RedBlackEntry class must define the member
// function GetKey().  The print() member function does not have to
// be defined because a default definition exists.
//
// The GetKey() function should return an integer key for that entry.
// The key for an entry should never change otherwise bad things might occur.

#define COLOR(color) (color == 1) ? "red" : "black"

static const short BLACK_ = 0;
static const short RED_ = 1; //2013.09.13 color is reversed compared with RBTreeC.h
static const short LEFT_ = 100;
static const short RIGHT_ = 200;

template<typename keyType, typename dataType>
class RedBlackTreeNode {
	//friend class RedBlackTree;
//protected:
	//RedBlackEntry * storedEntry;

public:
	keyType key;
	dataType* dataPtr;
	unsigned short color;
	//int key;
	//int red; /* if red=0 then the node is black */
	RedBlackTreeNode<keyType, dataType> * left;
	RedBlackTreeNode<keyType, dataType> * right;
	RedBlackTreeNode<keyType, dataType> * parent;

	RedBlackTreeNode() {
		parent = NULL;
		dataPtr = NULL;
		this->left = NULL;
		this->right = NULL;
		this->color = RED_;
	}
	RedBlackTreeNode(RedBlackTreeNode<keyType, dataType>* _parent, keyType _key,
			dataType* _dataPtr) :
			parent(_parent), key(_key), dataPtr(_dataPtr) {
		/*
		 * key_, data_ are references, have to be initialized in the away above.
		 */
		this->left = NULL;
		this->right = NULL;
		this->color = RED_;
	}
	//RedBlackTreeNode(RedBlackEntry *);
	//RedBlackEntry * GetEntry() const;
	~RedBlackTreeNode() {

	}
	RedBlackTreeNode<keyType, dataType> *getParent() {
		return this->parent;
	}
	RedBlackTreeNode<keyType, dataType> *getLeft() {
		return this->left;
	}
	RedBlackTreeNode<keyType, dataType> *getRight() {
		return this->right;
	}
	dataType* getDataPtr() {
		return this->dataPtr;
	}
	keyType getKey() {
		return this->key;
	}
	short getColor() {
		return this->color;
	}
	void setLeft(RedBlackTreeNode<keyType, dataType> *nodePtr) {
		this->left = nodePtr;
	}
	void setRight(RedBlackTreeNode<keyType, dataType> *nodePtr) {
		this->right = nodePtr;
	}
	void setParent(RedBlackTreeNode<keyType, dataType> *nodePtr) {
		this->parent = nodePtr;
	}
	void setDataPtr(dataType* dataPtr) {
		this->dataPtr = dataPtr;
	}
	void setKey(keyType key) {
		this->key = key;
	}

	void setColor(short color) {
		this->color = color;

	}

};

template<typename keyType, typename dataType>
class RedBlackTree {
protected:
	long _noOfNodes;	//size() is expensive, use this to keep track.
	short isValidRedBlackTreeRecur_(RedBlackTreeNode<keyType, dataType> *nodePtr) {
		if (nodePtr == nil)
			return 1;
		if (!isValidRedBlackTreeRecur_(nodePtr->getLeft())) {
			return 0;
		}
		if (!isValidRedBlackTreeRecur_(nodePtr->getRight())) {
			return 0;
		}
		if (nodePtr->getParent() == nil && nodePtr->getColor() == BLACK_) {
			return 1;
		}else if (nodePtr->getColor() == RED_){
			if (nodePtr->getParent() != nil){
				if (nodePtr->getParent()->getColor() == BLACK_) {
					return 1;
				}else{
					return 0;
				}
			}else if (nodePtr->getParent()==nil){
				return 1;
			}
		}else if (nodePtr->getColor() == BLACK_ && nodePtr->getParent() != nil
				&& nodePtr->getParent()->getColor() == BLACK_) {
			return 1;
		}
		return 0;
	}

	RedBlackTreeNode<keyType, dataType> *sibling_(
			RedBlackTreeNode<keyType, dataType> *nodePtr) {
		/*
		 * 2013.09.21 YH
		 */
		if (nodePtr->getParent() != nil) {
			if (nodePtr == nodePtr->getParent()->getLeft())
				return nodePtr->getParent()->getRight();
			else
				return nodePtr->getParent()->getLeft();
		} else {
			return nil;
		}
	}

	void leftRotate(RedBlackTreeNode<keyType, dataType>* x) {
		/***********************************************************************/
		/*  FUNCTION:  leftRotate */
		/**/
		/*  INPUTS:  the node to rotate on */
		/**/
		/*  OUTPUT:  None */
		/**/
		/*  Modifies Input: this, x */
		/**/
		/*  EFFECTS:  Rotates as described in _Introduction_To_Algorithms by */
		/*            Cormen, Leiserson, Rivest (Chapter 14).  Basically this */
		/*            makes the parent of x be to the left of x, x the parent of */
		/*            its parent before the rotation and fixes other pointers */
		/*            accordingly.  */
		/***********************************************************************/

		RedBlackTreeNode<keyType, dataType> * y;

		/*  I originally wrote this function to use the sentinel for */
		/*  nil to avoid checking for nil.  However this introduces a */
		/*  very subtle bug because sometimes this function modifies */
		/*  the parent pointer of nil.  This can be a problem if a */
		/*  function which calls leftRotate also uses the nil sentinel */
		/*  and expects the nil sentinel's parent pointer to be unchanged */
		/*  after calling this function.  For example, when DeleteFixUP */
		/*  calls leftRotate it expects the parent pointer of nil to be */
		/*  unchanged. */

		y = x->right;
		x->right = y->left;

		if (y->left != nil)
			y->left->parent = x; /* used to use sentinel here */
		/* and do an unconditional assignment instead of testing for nil */

		y->parent = x->parent;

		/* instead of checking if x->parent is the root as in the book, we */
		/* count on the root sentinel to implicitly take care of this case */
		if (x == x->parent->left) {
			x->parent->left = y;
		} else {
			x->parent->right = y;
		}
		y->left = x;
		x->parent = y;

#ifdef CHECK_RB_TREE_ASSUMPTIONS
		checkAssumptions();
#elif defined(DEBUG_ASSERT)
		Assert(!nil->color,"nil not red in RedBlackTree::leftRotate");
#endif
	}
	void rightRotate(RedBlackTreeNode<keyType, dataType>* y) {

		/***********************************************************************/
		/*  FUNCTION:  RighttRotate */
		/**/
		/*  INPUTS:  node to rotate on */
		/**/
		/*  OUTPUT:  None */
		/**/
		/*  Modifies Input?: this, y */
		/**/
		/*  EFFECTS:  Rotates as described in _Introduction_To_Algorithms by */
		/*            Cormen, Leiserson, Rivest (Chapter 14).  Basically this */
		/*            makes the parent of x be to the left of x, x the parent of */
		/*            its parent before the rotation and fixes other pointers */
		/*            accordingly.  */
		/***********************************************************************/

		RedBlackTreeNode<keyType, dataType> * x;

		/*  I originally wrote this function to use the sentinel for */
		/*  nil to avoid checking for nil.  However this introduces a */
		/*  very subtle bug because sometimes this function modifies */
		/*  the parent pointer of nil.  This can be a problem if a */
		/*  function which calls leftRotate also uses the nil sentinel */
		/*  and expects the nil sentinel's parent pointer to be unchanged */
		/*  after calling this function.  For example, when DeleteFixUP */
		/*  calls leftRotate it expects the parent pointer of nil to be */
		/*  unchanged. */

		x = y->left;
		y->left = x->right;

		if (nil != x->right)
			x->right->parent = y; /*used to use sentinel here */
		/* and do an unconditional assignment instead of testing for nil */

		/* instead of checking if x->parent is the root as in the book, we */
		/* count on the root sentinel to implicitly take care of this case */
		x->parent = y->parent;
		if (y == y->parent->left) {
			y->parent->left = x;
		} else {
			y->parent->right = x;
		}
		x->right = y;
		y->parent = x;

#ifdef CHECK_RB_TREE_ASSUMPTIONS
		checkAssumptions();
#elif defined(DEBUG_ASSERT)
		Assert(!nil->color,"nil not red in RedBlackTree::rightRotate");
#endif
	}
	void insertFix_(RedBlackTreeNode<keyType, dataType> *z) {
		/*  This function should only be called by RedBlackTree::Insert */

		/***********************************************************************/
		/*  FUNCTION:  insertFix_  */
		/**/
		/*  INPUTS:  z is the node to insert */
		/**/
		/*  OUTPUT:  none */
		/**/
		/*  Modifies Input:  this, z */
		/**/
		/*  EFFECTS:  Inserts z into the tree as if it were a regular binary tree */
		/*            using the algorithm described in _Introduction_To_Algorithms_ */
		/*            by Cormen et al.  This funciton is only intended to be called */
		/*            by the Insert function and not by the user */
		/***********************************************************************/

		RedBlackTreeNode<keyType, dataType> * x;
		RedBlackTreeNode<keyType, dataType> * y;

		z->left = z->right = nil;
		y = root;	//2013.09.22 YH now root is a regular node itself, rather than a proxy super node
		x = root->left;	//this is actual root
		while (x != nil) {
			y = x;
			if (x->key > z->key) {
				x = x->left;
			} else { /* x->key <= z->key */
				x = x->right;
			}
		}
		z->parent = y;
		/*
		if (y == nil){	//2013.09.22 YH, root is no longer a super node.
			root = z;
		} else
		*/
		if ((y == root) || (y->key > z->key)) {
			y->left = z;
		} else {
			y->right = z;
		}

#if defined(DEBUG_ASSERT)
		Assert(!nil->color,"nil not red in RedBlackTree::insertFix_");
#endif
	}

	void TreePrintHelper(RedBlackTreeNode<keyType, dataType> * x) const {
		if (x != nil) {
			TreePrintHelper(x->left);
			std::cout << x; //->print(nil, root);
			TreePrintHelper(x->right);
		}
	}
	void FixUpMaxHigh(RedBlackTreeNode<keyType, dataType> *) {
	}
	void deleteFixup_(RedBlackTreeNode<keyType, dataType> *x) {
		/***********************************************************************/
		/*  FUNCTION:  deleteFixup_ */
		/**/
		/*    INPUTS:  x is the child of the spliced */
		/*             out node in deleteNode. */
		/**/
		/*    OUTPUT:  none */
		/**/
		/*    EFFECT:  Performs rotations and changes colors to restore red-black */
		/*             properties after a node is deleted */
		/**/
		/*    Modifies Input: this, x */
		/**/
		/*    The algorithm from this function is from _Introduction_To_Algorithms_ */
		/***********************************************************************/

		RedBlackTreeNode<keyType, dataType> * w;
		RedBlackTreeNode<keyType, dataType> * rootLeft = root->left;

		while ((!x->color) && (rootLeft != x)) {
			if (x == x->parent->left) {
				w = x->parent->right;
				if (w->color) {
					w->color = BLACK_;
					x->parent->color = RED_;
					leftRotate(x->parent);
					w = x->parent->right;
				}
				if ((!w->right->color) && (!w->left->color)) {
					w->color = RED_;
					x = x->parent;
				} else {
					if (!w->right->color) {
						w->left->color = BLACK_;
						w->color = RED_;
						rightRotate(w);
						w = x->parent->right;
					}
					w->color = x->parent->color;
					x->parent->color = BLACK_;
					w->right->color = BLACK_;
					leftRotate(x->parent);
					x = rootLeft; /* this is to exit while loop */
				}
			} else { /* the code below is has left and right switched from above */
				w = x->parent->left;
				if (w->color) {
					w->color = BLACK_;
					x->parent->color = RED_;
					rightRotate(x->parent);
					w = x->parent->left;
				}
				if ((!w->right->color) && (!w->left->color)) {
					w->color = RED_;
					x = x->parent;
				} else {
					if (!w->left->color) {
						w->right->color = BLACK_;
						w->color = RED_;
						leftRotate(w);
						w = x->parent->left;
					}
					w->color = x->parent->color;
					x->parent->color = BLACK_;
					w->left->color = BLACK_;
					rightRotate(x->parent);
					x = rootLeft; /* this is to exit while loop */
				}
			}
		}
		x->color = BLACK_;

#ifdef CHECK_RB_TREE_ASSUMPTIONS
		checkAssumptions();
#elif defined(DEBUG_ASSERT)
		Assert(!nil->color,"nil not black in RedBlackTree::deleteFixup_");
#endif
	}
public:
	/*  A sentinel is used for root and for nil.  These sentinels are */
	/*  created when RedBlackTreeCreate is called.  root->left should always */
	/*  point to the node which is the root of the tree.  nil points to a */
	/*  node which should always be black but has arbitrary children and */
	/*  parent and no key or info.  The point of using these sentinels is so */
	/*  that the root and nil nodes do not require special cases in the code */
	RedBlackTreeNode<keyType, dataType> * root;
	RedBlackTreeNode<keyType, dataType> * nil;
	RedBlackTree() {
		nil = new RedBlackTreeNode<keyType, dataType>();
		nil->left = nil->right = nil->parent = nil;
		nil->color = BLACK_;
		//nil->key = MIN_INT;
		//nil->storedEntry = NULL;

		root = new RedBlackTreeNode<keyType, dataType>();
		root->parent = root->left = root->right = nil;
		root->color = BLACK_;
		//root->key = MAX_INT;
		//root->storedEntry = NULL;
		_noOfNodes = 0;
	}
	~RedBlackTree() {
		/*
		 RedBlackTreeNode<keyType, dataType> * x = root->left;
		 TemplateStack<RedBlackTreeNode<keyType, dataType> *> stuffToFree;

		 if (x != nil) {
		 if (x->left != nil) {
		 stuffToFree.Push(x->left);
		 }
		 if (x->right != nil) {
		 stuffToFree.Push(x->right);
		 }
		 // delete x->storedEntry;
		 delete x;
		 while (stuffToFree.NotEmpty()) {
		 x = stuffToFree.Pop();
		 if (x->left != nil) {
		 stuffToFree.Push(x->left);
		 }
		 if (x->right != nil) {
		 stuffToFree.Push(x->right);
		 }
		 // delete x->storedEntry;
		 delete x;
		 }
		 }
		 delete nil;
		 delete root;
		 */
	}
	void print() {
		TreePrintHelper(root->left);
	}
	int deleteNode(RedBlackTreeNode<keyType, dataType>* z) {
		/*  FUNCTION:  deleteNode */
		/**/
		/*    INPUTS:  tree is the tree to delete node z from */
		/**/
		/*    OUTPUT:  returns the RedBlackEntry stored at deleted node */
		/**/
		/*    EFFECT:  Deletes z from tree and but don't call destructor */
		/**/
		/*    Modifies Input:  z */
		/**/
		/*    The algorithm from this function is from _Introduction_To_Algorithms_ */
		/***********************************************************************/

		RedBlackTreeNode<keyType, dataType>* y;
		RedBlackTreeNode<keyType, dataType>* x;

		y = ((z->left == nil) || (z->right == nil)) ? z : getSuccessor_(z);
		x = (y->left == nil) ? y->right : y->left;
		if (root == (x->parent = y->parent)) { /* assignment of y->p to x->p is intentional */
			root->left = x;
		} else {
			if (y == y->parent->left) {
				y->parent->left = x;
			} else {
				y->parent->right = x;
			}
		}
		if (y != z) { /* y should not be nil in this case */

#ifdef DEBUG_ASSERT
			Assert( (y!=nil),"y is nil in deleteNode \n");
#endif
			/* y is the node to splice out and x is its child */

			y->left = z->left;
			y->right = z->right;
			y->parent = z->parent;
			z->left->parent = z->right->parent = y;
			if (z == z->parent->left) {
				z->parent->left = y;
			} else {
				z->parent->right = y;
			}
			if (!(y->color)) {
				y->color = z->color;
				deleteFixup_(x);
			} else
				y->color = z->color;
			delete z;
#ifdef CHECK_RB_TREE_ASSUMPTIONS
			checkAssumptions();
#elif defined(DEBUG_ASSERT)
			Assert(!nil->color,"nil not black in RedBlackTree::Delete");
#endif
		} else {
			if (!(y->color))
				deleteFixup_(x);
			delete y;
#ifdef CHECK_RB_TREE_ASSUMPTIONS
			checkAssumptions();
#elif defined(DEBUG_ASSERT)
			Assert(!nil->color,"nil not black in RedBlackTree::Delete");
#endif
		}
		_noOfNodes--;
		return 0;
	}
	RedBlackTreeNode<keyType, dataType> * insertNode(keyType key,
			dataType* dataPtr) {

		/*  Before calling InsertNode  the node x should have its key set */

		/***********************************************************************/
		/*  FUNCTION:  InsertNode */
		/**/
		/*  INPUTS:  newEntry is the entry to insert*/
		/**/
		/*  OUTPUT:  This function returns a pointer to the newly inserted node */
		/*           which is guarunteed to be valid until this node is deleted. */
		/*           What this means is if another data structure stores this */
		/*           pointer then the tree does not need to be searched when this */
		/*           is to be deleted. */
		/**/
		/*  Modifies Input: tree */
		/**/
		/*  EFFECTS:  Creates a node node which contains the appropriate key and */
		/*            info pointers and inserts it into the tree. */
		/***********************************************************************/

		RedBlackTreeNode<keyType, dataType> * y;
		RedBlackTreeNode<keyType, dataType> * x;
		RedBlackTreeNode<keyType, dataType> * newNode;

		x = new RedBlackTreeNode<keyType, dataType>(nil, key, dataPtr);
		/*
		if (root == nil) {
			x->left = x->right = nil;
			root = x;
			return x;
		} else {
		*/
		insertFix_(x);
		newNode = x;
		x->color = RED_;
		while (x->parent->color) { /* use sentinel instead of checking for root */
			if (x->parent == x->parent->parent->left) {
				y = x->parent->parent->right;
				if (y->color) {
					x->parent->color = BLACK_;
					y->color = BLACK_;
					x->parent->parent->color = RED_;
					x = x->parent->parent;
				} else {
					if (x == x->parent->right) {
						x = x->parent;
						leftRotate(x);
					}
					x->parent->color = BLACK_;
					x->parent->parent->color = RED_;
					rightRotate(x->parent->parent);
				}
			} else { /* case for x->parent == x->parent->parent->right */
				/* this part is just like the section above with */
				/* left and right interchanged */
				y = x->parent->parent->left;
				if (y->color) {
					x->parent->color = BLACK_;
					y->color = BLACK_;
					x->parent->parent->color = RED_;
					x = x->parent->parent;
				} else {
					if (x == x->parent->left) {
						x = x->parent;
						rightRotate(x);
					}
					x->parent->color = BLACK_;
					x->parent->parent->color = RED_;
					leftRotate(x->parent->parent);
				}
			}
		}
		root->left->color = BLACK_;	//2013.09.22 YH root->left is the real root. root is a sentinel
#ifdef CHECK_RB_TREE_ASSUMPTIONS
		checkAssumptions();
#elif defined(DEBUG_ASSERT)
		Assert(!nil->color,"nil not red in RedBlackTree::Insert");
		Assert(!root->color,"root not red in RedBlackTree::Insert");
#endif
		_noOfNodes++;
		return (newNode);
	}
	RedBlackTreeNode<keyType, dataType> * getPredecessor_(
			RedBlackTreeNode<keyType, dataType> * x) const {

		/***********************************************************************/
		/*  FUNCTION:  getPredecessor_  */
		/**/
		/*    INPUTS:  x is the node to get predecessor of */
		/**/
		/*    OUTPUT:  This function returns the predecessor of x or NULL if no */
		/*             predecessor exists. */
		/**/
		/*    Modifies Input: none */
		/**/
		/*    Note:  uses the algorithm in _Introduction_To_Algorithms_ */
		/***********************************************************************/
		RedBlackTreeNode<keyType, dataType> * y;

		if (nil != (y = x->left)) { /* assignment to y is intentional */
			while (y->right != nil) { /* returns the maximum of the left subtree of x */
				y = y->right;
			}
			return (y);
		} else {
			y = x->parent;
			while (x == y->left) {
				if (y == root)
					return (nil);
				x = y;
				y = y->parent;
			}
			return (y);
		}
	}
	RedBlackTreeNode<keyType, dataType> * getSuccessor_(
			RedBlackTreeNode<keyType, dataType> * x) const {
		/***********************************************************************/
		/*  FUNCTION:  getSuccessor_  */
		/**/
		/*    INPUTS:  x is the node we want the succesor of */
		/**/
		/*    OUTPUT:  This function returns the successor of x or NULL if no */
		/*             successor exists. */
		/**/
		/*    Modifies Input: none */
		/**/
		/*    Note:  uses the algorithm in _Introduction_To_Algorithms_ */
		/***********************************************************************/

		RedBlackTreeNode<keyType, dataType> * y;

		if (nil != (y = x->right)) { /* assignment to y is intentional */
			while (y->left != nil) { /* returns the minium of the right subtree of x */
				y = y->left;
			}
			return (y);
		} else {
			y = x->parent;
			while (x == y->right) { /* sentinel used instead of checking for nil */
				x = y;
				y = y->parent;
			}
			if (y == root)
				return (nil);
			return (y);
		}
	}
	RedBlackTreeNode<keyType, dataType> *queryTreeRecur_(
			RedBlackTreeNode<keyType, dataType> *nodePtr, keyType key) {
		if (nodePtr == nil) {
			return nil;
		}
		if (key < nodePtr->key) {
			return queryTreeRecur_(nodePtr->getLeft(), key);
		} else if (key > nodePtr->key) {
			return queryTreeRecur_(nodePtr->getRight(), key);
		} else {
			return nodePtr;
		}
	}
	RedBlackTreeNode<keyType, dataType> * queryTree(keyType key) {
		return queryTreeRecur_(root->left, key);
	}
	long size() {
		/*
		 * 2013.09.23 very expensive operation, use noOfNodes() to check.
		 */
		return count_(root->left, 0);
	}
	long noOfNodes(){
		return _noOfNodes;
	}
	long maxDepth() {
		return maxDepthRecur_(root->left);
	}
	RedBlackTreeNode<keyType, dataType> *getMinimum() {
		return getMinimum_(root->left);

	}
	RedBlackTreeNode<keyType, dataType> *getMaximum() {
		return getMaximum_(root->left);

	}

	long count_(RedBlackTreeNode<keyType, dataType> *nodePtr, long num) {
		if (nodePtr == nil) {
			return num;
		}
		return count_(nodePtr->getLeft(), count_(nodePtr->getRight(), ++num));
	}
	long maxDepthRecur_(RedBlackTreeNode<keyType, dataType> *nodePtr) {
		if (nodePtr == nil) {
			return 0;
		}
		long leftDepth = maxDepthRecur_(nodePtr->getLeft());
		long rightDepth = maxDepthRecur_(nodePtr->getRight());
		if (leftDepth > rightDepth) {
			return leftDepth + 1;
		} else {
			return rightDepth + 1;
		}
	}
	RedBlackTreeNode<keyType, dataType> *getMinimum_(
			RedBlackTreeNode<keyType, dataType> *nodePtr) {
		if (nodePtr == nil) {
			return nil;
		}
		if (nodePtr->getLeft() != nil) {
			return getMinimum_(nodePtr->getLeft());
		}
		return nodePtr;
	}
	RedBlackTreeNode<keyType, dataType> *getMaximum_(
			RedBlackTreeNode<keyType, dataType> *nodePtr) {
		if (nodePtr == nil) {
			return nil;
		}
		if (nodePtr->getRight() != nil) {
			return getMaximum_(nodePtr->getRight());
		}
		return nodePtr;
	}
	RedBlackTreeNode<keyType, dataType> *grandparent_(
			RedBlackTreeNode<keyType, dataType> *nodePtr) {
		if ((nodePtr != nil) && (nodePtr->getParent() != nil)) {
			return nodePtr->getParent()->getParent();
		} else {
			return nil;
		}
	}
	RedBlackTreeNode<keyType, dataType> *uncle_(
			RedBlackTreeNode<keyType, dataType> *nodePtr) {
		RedBlackTreeNode<keyType, dataType> *myGPNode = grandparent_(nodePtr);
		if (myGPNode == nil) {
			return nil;
		}
		if (nodePtr->getParent() == myGPNode->getLeft()) {
			return myGPNode->getRight();
		} else {
			return myGPNode->getLeft();
		}
	}
	//TemplateStack<RedBlackTreeNode<keyType, dataType> *> * Enumerate(int low, int high) ;
	void checkAssumptions() const {
		VERIFY(nil->color == 0);
		VERIFY(root->color == 0);
	}

	short isValidRedBlackTree() {
		return isValidRedBlackTreeRecur_(root->left);
	}

	bool isNULLNode(RedBlackTreeNode<keyType, dataType> * nodePtr){
		/*
		 * 2013.09.22 check whether a node is a leaf (nil)
		 */
		return nodePtr==nil;
	}

};

#endif
