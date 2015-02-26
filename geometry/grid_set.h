/***************************************************************************
 *            grid_set.h
 *
 *  Copyright  2008-12  Ivan S. Zapreev, Pieter Collins
 *
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Templece Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/*! \file grid_set.h
 *  \brief Grid paving is used to represent sets, based on integer and dyadic coordinate cells, of a grid.
 */

#ifndef ARIADNE_GRID_SET_H
#define ARIADNE_GRID_SET_H

#include <iostream>
#include <iomanip>
#include <string>

#include <boost/iterator/iterator_facade.hpp>
#include <memory>

#include "utility/tribool.h"
#include "utility/array.h"

#include "utility/binary_word.h"

#include "utility/exceptions.h"
#include "geometry/box.h"
#include "geometry/point.h"
#include "geometry/list_set.h"

#include "numeric/numeric.h"

#include "geometry/set_interface.h"
#include "geometry/paving_interface.h"
#include "algebra/vector.h"
#include "geometry/grid.h"
#include "geometry/grid_cell.h"

#include "output/graphics_interface.h"

using namespace std;
using namespace Ariadne;

namespace Ariadne {

/*Some type definitions*/
typedef std::vector<Bool> BooleanArray;
typedef Array<Int> IndexArray;
typedef Array<SizeType> SizeArray;

typedef unsigned short dimension_type;

/*Some pre-declarations*/
class BinaryTreeNode;
class Grid;
class GridAbstractCell;
class GridCell;
class GridOpenCell;
class GridTreeSubset;
class GridTreeSet;

class GridTreeCursor;
class GridTreeConstIterator;

/*Declarations of classes in other files*/
template<class BS> class ListSet;

OutputStream& operator<<(OutputStream& output_stream, const BinaryTreeNode & binary_tree );
OutputStream& operator<<(OutputStream& os, const GridCell& theGridCell);
OutputStream& operator<<(OutputStream& os, const GridOpenCell& theGridOpenCell );
OutputStream& operator<<(OutputStream& os, const GridTreeCursor& theGridTreeCursor);
OutputStream& operator<<(OutputStream& os, const GridTreeSubset& theGridTreeSubset);
OutputStream& operator<<(OutputStream& os, const GridTreeSet& theGridTreeSet);

Bool subset(const GridCell& theCell, const GridTreeSubset& theSet);
Bool intersect(const GridCell& theCell, const GridTreeSubset& theSet);
Bool subset(const GridTreeSubset& theSet1, const GridTreeSubset& theSet2);
Bool intersect(const GridTreeSubset& theSet1, const GridTreeSubset& theSet2);

GridTreeSet join(const GridTreeSubset& theSet1, const GridTreeSubset& theSet2);
GridTreeSet intersection(const GridTreeSubset& theSet1, const GridTreeSubset& theSet2);
GridTreeSet difference(const GridTreeSubset& theSet1, const GridTreeSubset& theSet2);

GridTreeSet outer_approximation(const ExactBox& theBox, const Grid& theGrid, const Nat numSubdivInDim);
GridTreeSet outer_approximation(const CompactSetInterface& theSet, const Grid& theGrid, const Nat numSubdivInDim);
GridTreeSet outer_approximation(const CompactSetInterface& theSet, const Nat numSubdivInDim);
template<class BS> GridTreeSet outer_approximation(const ListSet<BS>& theSet, const Nat numSubdivInDim);
GridTreeSet inner_approximation(const OpenSetInterface& theSet, const Grid& theGrid, const Nat height, const Nat numSubdivInDim);

#ifdef ARIADNE_ENABLE_SERIALIZATION
  template<class A> Void serialize(A& archive, const GridTreeSet& set, const unsigned int version);
#endif /* ARIADNE_ENABLE_SERIALIZATION */


/*! \brief The binary-tree node operation is not allowed on a non-leaf node. */
class NotALeafNodeException : public std::logic_error {
  public:
    NotALeafNodeException(const StringType& str) : std::logic_error(str) { }
};

/*! \brief The binary-tree node operation is not allowed on a leaf node. */
class IsALeafNodeException : public std::logic_error {
  public:
    IsALeafNodeException(const StringType& str) : std::logic_error(str) { }
};

/*! \brief The GridTreeCursor throws this exception if we try to go beyond the binary tree. */
class NotAllowedMoveException : public std::logic_error {
  public:
    NotAllowedMoveException(const StringType& str) : std::logic_error(str) { }
};


/*! \brief The binary tree node.
 *
 * This node is to be used in a binary tree designed for subdividing the state
 * space into enabled and disabled cells. This is required for representing
 * subsets of the state space.
 *
 * \b Storage: We only store pointers to the left and right subtrees and the Kleenean
 * value indicating whether this cell is enabled/disabled or we do not know.
 */
class BinaryTreeNode {
  protected:
    /*! \brief Defines whether the given node of the tree is on/off or we do not know*/
    Kleenean _isEnabled;

    /*! \brief The left and right subnodes of the tree. Note that,
     * \a pLeftNode == \a NULL iff \a pRightNode == \a NULL. The latter
     * is allowed iff \a _isEnabled == \a TRUE or \a _isEnabled == \a FALSE,
     * i.e. the node can be a leaf.
     */
    BinaryTreeNode* _pLeftNode;
    BinaryTreeNode* _pRightNode;

    /*! \brief This method splits the enabled subtrees of the tree rooted to
     * \a pCurrentNode in such a way that the tree depth becomes \a depth.
     * If the initial tree depth is greater than \a depth then nothing is done.
     */
    Void mince_node(BinaryTreeNode* pCurrentNode, const Nat depth);

    /*! \brief This method recombined the sub tree nodes rooted to \a pCurrentNode.
     * Note that, the two leaf nodes with the same parent are removed if they have
     * the same value of isEnabled fields.
     */
    Void recombine_node(BinaryTreeNode * pCurrentNode);

    /*! \brief This method is used for recursive restoration of the binary tree from
     *  the used data, i.e. \a theTree and \a theEnabledCells. It is used in the
     *  constructor \a BinaryTreeNode( const BooleanArray& , const BooleanArray& )
     */
    Void restore_node( BinaryTreeNode * pCurrentNode, Nat & arr_index, Nat & leaf_counter,
                       const BooleanArray& theTree, const BooleanArray& theEnabledCells);

    /*! \brief This method is used in constructors for the node initialization */
    Void init( Kleenean isEnabled, BinaryTreeNode* pLeftNode, BinaryTreeNode* pRightNode );

  public:
    //@{
    //! \name Constructors

    /*! \brief Construct a tree node. */
    explicit BinaryTreeNode(const Kleenean _isEnabled = false );

    /*! \brief The copy constructor.
     * The the node and all it's sub nodes are copied.
     */
    explicit BinaryTreeNode(const BinaryTreeNode& theTreeNode);

    /*! \brief Constructs a binary tree from the boolean arrays. \a theTree defines
     * the tree structure, \a theEnabledCells defines the tree nodes.
     * IVAN S. ZAPREEV:
     * WARNING: We assume that every node has either no children or both of them!
     * NOTE: The dinary data of the tree is organized in the following way:
     *          1      The \a theTree contains Depth first search lay out of the tree,
     *         / \     where 1 stands for the non-leaf node and 0 for a leaf node, we
     *        /   \    always visit the left su-node first, e.g. the tree on the left
     *       2     5   is encodes the Array:     [1, 1, 0, 0, 0]
     *      / \        where the corresponding    ^  ^  ^  ^  ^
     *     /   \       tree nodes are:            1  2  3  4  5
     *    3     4      \a theEnabledCells contains true/false values for the leaf nodes
     * of the tree. Their order is the same as in \a theTree, e.g. here it is: 3,4,5
     */
    explicit BinaryTreeNode( const BooleanArray& theTree, const BooleanArray& theEnabledCells );

    //@}

    ~BinaryTreeNode();

    //@{
    //! \name Properties

    /*! \brief Returns true if the node is marked as enabled, otherwise false */
    Bool is_enabled() const;

    /*! \brief Returns true if some of the leaf nodes in the tree rooted to this node are enabled, otherwise false */
    Bool has_enabled() const;

    /*! \brief Returns true if all leaf nodes in the tree rooted to this node are enabled, otherwise false */
    Bool all_enabled() const;

    /*! \brief This method returns true if the given path defines a node in the tree that is either enabled
     *  or is in a "virtual" subtree of some enabled node (note that enabled nodes can only be leafs).
     * Note that: \a path is treated as if it is rooted to this node, we assume that the path starts
     * from position \a position of \a path. The parameter \a position is used for recursive calls only.
     */
    Bool is_enabled( const BinaryWord & path, const Nat position = 0) const;

    /*! \brief Returns true if the node is marked as disabled, otherwise false */
    Bool is_disabled() const;

    /*! \brief Returns true if the node is a leaf (pLeftNode == NULL && pRightNode == NULL) otherwise false */
    Bool is_leaf() const;

    /*! \brief Return the left or right sub-node */
    BinaryTreeNode * child_node(Bool left_or_right) const;

    /*! \brief Return the left sub-node */
    BinaryTreeNode * left_node() const;

    /*! \brief Return the right sub-node */
    BinaryTreeNode * right_node() const;

    /*! \brief Returns the depth of the sub-tree rooted to the given node, i.e. the depth of it's deepest node */
    Nat depth() const;

    /*! \brief Allows to compare to binaty tree nodes */
    Bool operator==(const BinaryTreeNode & otherNode ) const;

    /*! \brief Marks the node as enabled or disabled. */
    Void set(Bool all_enabled_or_all_disabled);

    static Bool is_equal_nodes( const BinaryTreeNode * pFirstNode, const BinaryTreeNode * pSecondNode );

    //@}

    //@{
    //! \name Leaf Operations

    /*! \brief This method makes the node to become a leaf node with the enabled value : \a is_enabled
     * NOTE: the leat and the right sub-trees (are deallocated.
     * WARNING: this method MUST NOT be called on a non-leaf node!!!
     */
    Void make_leaf( Kleenean is_enabled );

    /*! \brief Marks the leaf node as enabled, otherwise they through \a NotALeafNodeEsception */
    Void set_enabled();

    /*! \brief Marks the leaf node as disabled, otherwise they through \a NotALeafNodeEsception */
    Void set_disabled();

    /*! \brief Marks the node as neither enabled nor disabled, is only applicable to non-leaf nodes.
     * When applied to a leaf node, throws IsALeafNodeException.
     */
    Void set_unknown();

    /*! \brief Splits the leaf node, i.e. adds two subnodes with the _isEnabled field value inherited from the parent node.
     * The parent's _isEnabled is set to intermediate, because the subsequent operation on the subtree might enable/disable
     * subnodes and we do not want to keep track of these changes. If the node is not a leaf then nothing is done.
     */
    Void split();

    /*! \brief This method splits the enabled subtrees of the tree rooted to
     * this node in such a way that the tree depth becomes \a depth.
     * If the initial tree depth is greater than \a depth then nothing is done.
     */
    Void mince(const Nat depth);

    //@}

    //@{
    //! \name

    /*! \brief Allows to assign one binary tree node to the other, this is done by copying
     * the nodes value and the sub-trees. The sub-trees of this node are deallocated.
     */
    BinaryTreeNode& operator=(const BinaryTreeNode & otherNode );

    /*! \brief Copy all the data (including the sub-nodes) from the node pointed by \a pOtherNode into the given node.
     * Note that here we will create copies of the sub nodes, and NOT just copy pointers to them!
     */
    Void copy_from( const BinaryTreeNode * pOtherNode );

    /*! \brief This method recombined the sub tree nodes. Note that, the two leaf nodes with
     * the same parent are removed if they have the same value of isEnabled fields.
     */
    Void recombine();

    /*! \brief Stores the binary tree in a form of two arrays, their structure is the same as needed for
     *   the BinaryTreeNode( const BooleanArray& , const BooleanArray&  ) constructor
     */
    Void tree_to_binary_words( BinaryWord & tree, BinaryWord & leaves ) const;

    /*! \brief Stores the binary tree node as a string*/
    string node_to_string() const;

    /*! \brief Finds(creates) the leaf node defined by the \a path and marks it as enabled.
     * If some prefix of the \a path references an enabled node then nothing is done.
     */
    Void add_enabled( const BinaryWord & path );

    /*! \brief This method adjoins the enabled nodes of \a subTree to this tree.
     * Note that, the position of the root node of \a subTree within this
     * tree is defined by path.
     */
    Void add_enabled( const BinaryTreeNode * pOtherSubTree, const BinaryWord & path );

    /*! \brief This method merges \a pFromTreeRoot into \a pToTreeRoot.
     *  The enabled nodes of the former tree are added to the latter tree.
     *  If in the latter tree there is an enabled leaf node and in the former
     *  tree the corresponding node is not a leaf, then this node of the latter
     *  tree stays intact. If on the other hand we have a non-leaf node in
     *  \a pToTreeRoot and we are adding to it an enabled node of \a pFromTreeRoot
     *  then we just "substitute" the former one with the latter.
     *  NOTE: 1. This function is recursive. 2. No pointers are copied between
     *  \a pToTreeRoot and \a pFromTreeRoot.
     */
    Void add_enabled( BinaryTreeNode* pToTreeRoot, const BinaryTreeNode* pFromTreeRoot );

    /*! \brief Starting in the \a pNode node as at the root, this method counts
     *  the number of enabled leaf nodes in the subtree
     *  rooted at pNode.
     */
    static SizeType count_enabled_leaf_nodes( const BinaryTreeNode* pNode );

    /*! \brief Starting in the \a pRootTreeNode node as at the root, this method finds(creates)
     *  the leaf node defined by the \a path and marks it as enabled. If some prefix of the \a path
     *  references an enabled node then nothing is done.
     *  NOTE: This is a recursive method on the position (\a position) in the binary path (\a path).
     *  Therefore, the initial evaluate of this method should be done with \a position == 0;
     */
    static Void add_enabled( BinaryTreeNode* pRootTreeNode, const BinaryWord& path, const Nat position = 0 );

    /*! \brief Creates a binary tree of the height rootNodePath.size(), puts the subtree oldRootNode
     * into the node defined by the path \a rootNodePath, returns the root node of the extended tree.
     */
    static BinaryTreeNode * prepend_tree( const BinaryWord & rootNodePath, BinaryTreeNode * oldRootNode);

    /*! \brief This method restricts \a pThisNode to \a pOtherNode.
     * In essance we do the inplace AND on the tree node pThisNode.
     * Note that, this method is recursive.
     */
    static Void restrict( BinaryTreeNode * pThisNode, const BinaryTreeNode * pOtherNode );

    /*! \brief This method removed enabled nodes of \a pOtherNode from \a pThisNode.
     * Note that, this method is recursive.
     */
    static Void remove( BinaryTreeNode * pThisNode, const BinaryTreeNode * pOtherNode );

    /*! \brief checks if two trees intersect in a set-theory sence.
     * I.e. we assume that pRootNodeOne and pRootNodeTwo correspond to the same (virtual) root
     * and then we see if the enabled leaf node of one tree contain enabled leaf nodes of
     * another tree as their (virtual) children.
     */
    static Bool intersect( const BinaryTreeNode * pRootNodeOne, const BinaryTreeNode * pRootNodeTwo );

    /*! \brief checks if the tree pRootNodeOne is a subset of the tree pRootNodeTwo, in a set-theory sence.
     * I.e. we assume that pRootNodeOne and pRootNodeTwo correspond to the same (virtual) root
     * and then we see if every enabled leaf node of pRootNodeOne is contained in the enabled leaf nodes of
     * pRootNodeTwo, or it is covered by the enabled leaf nodes of pRootNodeTwo. When we write, contained and
     * covered then we mean: is a subnode in the (virtual) tree and all it's subnodes in the (virtual) tree.
     */
    static Bool subset( const BinaryTreeNode * pRootNodeOne, const BinaryTreeNode * pRootNodeTwo );

    //@}
};


/*! \brief This class represents a subpaving of a paving. Note that, the subtree enclosed into
 * this class is just a pointer to the node in the tree of some paving. This class is not
 * responsible for deallocation of that original tree.
 */
class GridTreeSubset
    : public virtual SubPavingInterface
    , public virtual DrawableInterface
{
  protected:

    friend class GridTreeCursor;

    friend GridTreeSet outer_approximation( const CompactSetInterface& theSet, const Grid& theGrid, const Nat numSubdivInDim );
    friend GridTreeSet inner_approximation( const OpenSetInterface& theSet, const Grid& theGrid, const Nat numSubdivInDim );

    /*! \brief The pointer to the root node of the subpaving tree.
     * Note that, this is not necessarily the root node of the corresponding paving tree.
     */
    BinaryTreeNode * _pRootTreeNode;

    /*! \brief The paving cell corresponding to the root node of the SubPaving.*/
    GridCell _theGridCell;

    /*! \brief this function takes the interval width and computes how many binary subdivisions
     * one has to make in order to have sub-intervals of the width <= \a theMaxWidth
     */
    Nat compute_number_subdiv( Float64 theWidth, const Float64 theMaxWidth) const;

    /*! \brief This method checks whether the set defined by \a pCurrentNode is a superset
     *  of \a theBox, in case when it is known that the cell corresponding to the root of
     *  pCurrentNode [and defined by \a theGrid, the primary cell (\a theHeight) to which
     *  this tree is (virtually) rooted via the path theWord] encloses \a theBox.
     *  This is a recursive procedure and it returns true only if there are no disabled
     *  cells in \a pCurrentNode that intersect with theBox.
     */
    static Kleenean superset( const BinaryTreeNode* pCurrentNode, const Grid& theGrid,
                             const Nat theHeight, BinaryWord &theWord, const ExactBox& theBox );

    /*! \brief This method checks whether the set defined by \a pCurrentNode is a subset
     *  of \a theBox. The set of \a pCurrentNode is defined by \a theGrid, the primary
     *  cell (\a theHeight) to which this tree is (virtually) rooted via the path theWord.
     *  This is a recursive procedure and it returns true only if all enabled sub-cells of
     *  \a pCurrentNode are sub-sets of \a theBox.
     */
    static Kleenean subset( const BinaryTreeNode* pCurrentNode, const Grid& theGrid,
                           const Nat theHeight, BinaryWord &theWord, const ExactBox& theBox );

    /*! \brief This method checks whether \a theBox is disjoint from the set defined by
     *  \a pCurrentNode, \a theGrid, the primary cell (\a theHeight) to which this
     *  tree is (virtually) rooted via the path theWord. This is done using the recursive
     *  procedure by checking the cells of the tree that intersect with the box and going
     *  down to the leaves. When reaching a leaf node that is enabled we conclude that we
     *  have an intersection. If there are no such nodes then there is no intersection,
     *  and the sets are disjoint.
     */
    static Kleenean disjoint( const BinaryTreeNode* pCurrentNode, const Grid& theGrid,
                             const Nat theHeight, BinaryWord &theWord, const ExactBox& theBox );

    /*! \brief This method checks whether \a theBox overlaps the set defined by
     *  \a pCurrentNode, \a theGrid, the primary cell (\a theHeight) to which this
     *  tree is (virtually) rooted via the path theWord. This is done using the recursive
     *  procedure by checking the cells of the tree that overlap the box and going
     *  down to the leaves. When reaching a leaf node that is enabled we conclude that we
     *  have an overlap. If there are no such nodes then there is no overlap. Note that this
     *  method cannot be implemented using disjoint since disjoint tests if the closures
     *  are disjoint, and overlaps tests if the interiors are not disjoint.
     */
    static Kleenean intersects( const BinaryTreeNode* pCurrentNode, const Grid& theGrid,
                               const Nat theHeight, BinaryWord &theWord, const ExactBox& theBox );

    /*! Allow to convert the number of subdivisions in each dimension, i.e. \a numSubdivInDim,
     *  starting from the zero cell into the number of subdivisions that have to be done in a
     *  subpaving rooted to the primary cell of the height \a primaryCellHeight and having a
     *  length of the path from the primary cell to the root cell of length \a primaryToRootCellPathLength.
     *  Returns zero if the binary tree already has a sufficient number of subdivisions due to
     *  its root cell being smaller than the cells obtained when subdividing \a numSubdivInDim
     *  times in each dimension starting from the zero-cell level.
     */
    Int zero_cell_subdivisions_to_tree_subdivisions( const Nat numSubdivInDim, const Nat primaryCellHeight,
                                                     const Nat primaryToRootCellPathLength ) const;

  public:

    /*! \brief A short name for the constant Iterator */
    typedef GridTreeConstIterator ConstIterator;

    //@{
    //! \name Constructors

    /*! \brief The new root node can only be constructed from the existing tree node.
     * Here, \a pRootTreeNode is not copied, we simply store its pointer inside this class.
     * Note that, \a pRootTreeNode should correspond to the sub-paving root node. Thus,
     * \a theHeight defines the height of the primary root cell of the GridTreeSet
     * (Remember that every GridTreeSubset is just a reference to a subtree of a GridTreeSet).
     * \a theWord defines the path to the \a pRootTreeNode node from the primary root cell
     * of the corresponding GridTreeSet.
     */
    GridTreeSubset( const Grid& theGrid, const Nat theHeight, const BinaryWord& theWord, BinaryTreeNode * pRootTreeNode );

    /*! \brief A copy constructor that only copies the pointer to the root of the binary tree and the cell */
    GridTreeSubset( const GridTreeSubset &otherSubset);

    /*! \brief Make a dynamically-allocated copy as a GridTreeSet. Required for DrawableInterface. */
    GridTreeSubset* clone() const;

    //@}

    /*! Virtual destructor. The destructor needs to be virtual since GridTreeSet is a subclass
     *  with different memory management. */
    virtual ~GridTreeSubset();

    //@{
    //! \name Properties

    /*! \brief True if the set is empty. */
    Bool is_empty() const;

    /*! \brief The number of activated cells in the set. */
    SizeType size() const;

    /*! \brief The dimension of the set. */
    DimensionType dimension() const;

    /*! \brief Returns a constant reference to the underlying grid. */
    const Grid& grid() const;

    /*! \brief Returns the const pointer to the root BinaryTreeNode of the SubPaving*/
    const BinaryTreeNode * binary_tree() const;

    /*! Recalculate the depth of the tree rooted at \a _pRootTreeNode */
    Nat depth() const;

    /*! The measure (area, volume) of the set in Euclidean space. */
    double measure() const;

    /*! The a branch along the binary tree. */
    GridTreeSubset branch(Bool left_or_right) const;

    /*! \brief Returns the \a GridCell corresponding to the ROOT NODE of this \a GridTreeSubset
     * WARNING: It is NOT the primary cell of the paving!
     */
    GridCell root_cell() const;

    /*! \brief Returns the \a GridCell corresponding to the ROOT NODE of \em this \a GridTreeSubset
     */
    GridCell cell() const { return this->root_cell(); };

    /*! \brief Computes a bounding box for a grid set. */
    UpperBox bounding_box() const;

    /*! \brief Allows to test if the two subpavings are "equal". The method returns true if
     * the grida are equal and the binary trees are equal. Note that, only in case both
     * GridTreeSubset objects are recombines, this method is guaranteed to tell you that
     * the two GridTreeSubset represent equal sets.
     */
    Bool operator==(const GridTreeSubset& anotherGridTreeSubset) const;

    //@}

    //@{
    //! \name Modifying operations

    /*! \brief Sets the ROOT NODE of this \a GridTreeSubset to either enabled (true) or disabled (false).
     */
    Void set_root_cell(Bool enabled_or_disabled);

    //@}

    //@{
    //! \name Subdivisions

    /*! \brief Subdivides the tree in such a way thaty it's depth becomes ( height + numSubdivInDim ) * D
     * Where height is the height of the primary cell to which the tree is rooted, \a numSubdivInDim is
     * the number of subdivisions in each dimension, and D is the number of dimensions of our space.
     * Note that, in case the subset is already subdivided to the required depth then nothing is done.
     * The latter can happen if the root cell of the subset is below the depth ( height + numSubdivInDim ) * D.
     */
    Void mince( Nat numSubdivInDim );

    /*! \brief Subdivides the tree up to the depth specified by the parameter.
     * Note that, we start from the root of the sub-paving and thus the subdivision
     * is done relative to it but not to the root of the original paving.
     */
    Void mince_to_tree_depth( const Nat theNewDepth );

    /*! \brief Subdivide the paving until the smallest depth such that the leaf
     * cells size is <= \a theMaxCellWidth. Note that, the disabled cells are
     * not subdivided.
     */
    Void subdivide( Float64 theMaxCellWidth );

    /*! \brief Recombines the subdivisions, for instance if all subcells of a cell are
     * enabled/disabled then they are put together.
     */
    Void recombine();

    //@}

    //@{
    //! \name Geometric Predicates

    /*! \brief Tests if a cell is a subset of the set. */
    Bool superset( const GridCell& theCell) const { return Ariadne::subset(theCell,*this); }

    /*! \brief Tests if a cell is a subset of a set. */
    friend Bool subset( const GridCell& theCell, const GridTreeSubset& theSet );

    /*! \brief Tests if a cell intersect (as an open set) a paving set. */
    friend Bool intersect( const GridCell& theCell, const GridTreeSubset& theSet );

    /*! \brief Tests if a grid set \a theSet1 is a subset of \a theSet2. */
    friend Bool subset( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 );

    /*! \brief Tests if two grid paving sets intersect (as open sets)
     *  If at least one of the GridTreeSubsets represents an empty set, then the result is false.
     */
    friend Bool intersect( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 );

    /*! \brief Tests if a grid set equals another paving. */
    virtual Bool equals(const SubPavingInterface&) const;

    /*! \brief Tests if a grid set is a subset of another paving. */
    virtual Bool subset(const SubPavingInterface&) const;

    /*! \brief Tests if a grid set intersects another paving. */
    virtual Bool intersects(const SubPavingInterface&) const;

    /*! \brief Tests if a grid set is a subset of a box. */
    Bool subset( const ExactBox& theBox ) const;

    /*! \brief Tests if a grid set is a superset of a box. */
    Bool superset( const ExactBox& theBox ) const;

    /*! \brief Tests if a grid set intersects (the interior) of a box. */
    Bool intersects( const ExactBox& theBox ) const;

    /*! \brief Tests if a grid set is disjoint from (the interior of) box. */
    Bool disjoint( const ExactBox& theBox ) const;

    /*! \brief Tests if the interior of a grid set is a superset of a box. */
    Sierpinski covers( const ExactBox& theBox ) const;

    /*! \brief Tests if (the closure of) a grid set is a subset of the interior of box. */
    Sierpinski inside( const ExactBox& theBox  ) const;

    /*! \brief Tests if (the closure of) a grid set is disjoint from (the closure of) a box. */
    Sierpinski separated( const ExactBox& theBox  ) const;

    /*! \brief Tests if a grid set overlaps (intersects the interior of) a box. */
    Sierpinski overlaps( const ExactBox& theBox ) const;

    //@}

    //@{
    //! \name Iterators

    /*! \brief A constant Iterator through the enabled leaf nodes of the subpaving. */
    ConstIterator begin() const;

    /*! \brief A constant Iterator to the end of the enabled leaf nodes of the subpaving. */
    ConstIterator end() const;

    //@}

    //@{
    //! \name Conversions

    /*! \brief An assignment operator that only copies the pointer to the root of the binary tree and the cell */
    GridTreeSubset& operator=( const GridTreeSubset &otherSubset);

    /*! \brief Convert to a list of ordinary boxes, unrelated to the grid. */
    operator ListSet<ExactBox>() const;

    //@}

    //@{
    //! \name Input/Output

    /*! \brief Draw on a two-dimensional canvas. */
    Void draw(CanvasInterface& canvas, const Projection2d& projection) const;

    /*! \brief Write to an output stream. */
    OutputStream& write(OutputStream& os) const;
    //@}

  private:
    virtual GridTreeSubset* _branch(Bool left_or_right) const;
    virtual ForwardConstantIteratorInterface<GridCell>* _begin() const;
    virtual ForwardConstantIteratorInterface<GridCell>* _end() const;

};

/*! \ingroup ListSetSubModule
 * \ingroup StorageModule
 * \brief The %GridTreeSet class represents a set of cells with mixed integer and dyadic coordinates.
 * The cells can be enabled or disabled (on/off), indicating whether they belong to the paving or not.
 * It is possible to have cells that are neither on nor off, indicating that they have enabled and
 * disabled sub cells.
 *
 * A trivial cell (level 0) of the paving corresponds to a unit hypercube, of the n-dimensional state
 * space. Subdivision of any 0-level paving is based on dyadic numbers and thus is a binary-style
 * partitioning of the cell. All cells with purely integer coordinates are enclosed in a (virtual)
 * binary tree. We illustrate this for a two dimensional case:
 *   1. Take the unit cell [0, 0]*[1, 1],
 *   2. Cells [-1, 0]*[0, 1] and [0, 0]*[1, 1] are rooted to [-1, 0]*[1, 1],
 *   3. Cells [-1, -1]*[0, 0] and [0, -1]*[1, 0] are rooted to [-1, -1]*[1, 0],
 *   4. Cells [-1, -1]*[1, 0] and [-1, 0]*[1, 1] are rooted to [-1, -1]*[1, 1].
 *
 * \sa GridCell
 */
class GridTreeSet
    : public virtual PavingInterface
    , public GridTreeSubset
{

  public:

    typedef GridCell CellType;

  protected:

    friend class GridTreeCursor;

    /*! \brief This method takes the height of the primary cell
     *  \a otherPavingPCellHeight and if it is:
     *    (a) higher then for this paving, pre-pends \a _pRootTreeNode.
     *    (b) lower then for this paving, locates it in this paging's tree.
     *    (c) equal to the hight of this paving's primary cell, does nothing.
     *  This method returns the node corresponding to the primary cell of
     *  height \a otherPavingPCellHeight in the (updated) paving.
     *  If \a stop_on_enabled is set to true then for the case (b) if we meet
     *  a leaf node on the path from the primary node of the paving to the primary
     *  node of the cell, we stop locating the node corresponding to the primary
     *  cell of height \a otherPavingPCellHeight and set \a has_stopped to true.
     *  If \a stop_on_disabled is set to true then for the case (b) if we meet
     *  a leaf node on the path from the primary node of the paving to the primary
     *  node of the cell, we stop locating the node corresponding to the primary
     *  cell of height \a otherPavingPCellHeight and set \a has_stopped to true.
     */
    BinaryTreeNode* align_with_cell( const Nat otherPavingPCellHeight, const Bool stop_on_enabled, const Bool stop_on_disabled, Bool & has_stopped );

    /*! \brief This method adjoins the outer approximation of \a theSet (computed on the fly) to this paving.
     *  We use the primary cell (enclosed in this paving) of height \a primary_cell_hight and represented
     *  by the paving's binary node \a pBinaryTreeNode. When adding the outer approximation, we compute it
     *  up to the level of accuracy given by \a max_mince_depth. This parameter defines, how many subdivisions
     *  of the binary tree we should make to get the proper cells for outer approximating \a theSet.
     *  This method is recursive, the parameter \a pPath defines the path to the current node pBinaryTreeNode
     *  from the root node in recursive calls, thus the initial evaluate for this method must be done with an empty word.
     */
    static Void _adjoin_outer_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const Nat primary_cell_height,
                                             const Nat max_mince_depth, const CompactSetInterface& theSet, BinaryWord * pPath );

    /*! \brief This method adjoins the inner approximation of \a theSet (computed on the fly) to this paving.
     *  We use the primary cell (enclosed in this paving) of height \a primary_cell_height and represented
     *  by the paving's binary node \a pBinaryTreeNode. When adding the inner approximation, we compute it
     *  up to the level of accuracy given by \a max_mince_depth. This parameter defines, how many subdivisions
     *  of the binary tree we should make to get the proper cells for inner approximating \a theSet.
     *  This method is recursive, the parameter \a pPath defines the path to the current node pBinaryTreeNode
     *  from the root node in recursive calls, thus the initial evaluate for this method must be done with an empty word.
     */
    static Void _adjoin_inner_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const Nat primary_cell_height,
                                             const Nat max_mince_depth, const OpenSetInterface& theSet, BinaryWord * pPath );

    /*! \brief This method adjoins the lower approximation of \a theSet (computed on the fly) to this paving.
     *  We use the primary cell (enclosed in this paving) of height \a primary_cell_hight and represented
     *  by the paving's binary node \a pBinaryTreeNode. When adding the lower approximation, we compute it
     *  up to the level of accuracy given by \a max_mince_depth. This parameter defines, how many subdivisions
     *  of the binary tree we should make to get the proper cells for lower approximating \a theSet.
     *  This method is recursive, the parameter \a pPath defines the path to the current node pBinaryTreeNode
     *  from the root node in recursive calls, thus the initial evaluate for this method must be done with an empty word.
     *  The approximation method does not recombine cells, as knowing that both children intersect a set is more
     *  information than knowing that the parent does.
     */
    static Void _adjoin_lower_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const Nat primary_cell_height,
                                             const Nat max_mince_depth, const OvertSetInterface& theSet, BinaryWord * pPath );

    /*! \brief This method adjoins the lower approximation of \a theSet (computed on the fly) to this paving.
     *  It is specialised for open sets, for which we have the superset() operator. If a set is a superset of
     *  a cell, then we know it overlaps the cell and all its children.
     */
    static Void _adjoin_lower_approximation( const Grid & theGrid, BinaryTreeNode * pBinaryTreeNode, const Nat primary_cell_height,
                                             const Nat max_mince_depth, const OpenSetInterface& theSet, BinaryWord * pPath );

    /*! \brief This method is uset to do restriction of this set to the set given by
     *  \a theOtherSubPaving Note that, here we require that the height of the primary
     *  root cell of this set is >= the height of the primary root cell of \a theOtherSubPaving.
     */
    Void restrict_to_lower( const GridTreeSubset& theOtherSubPaving );

    /*! \brief This method is uset to remove \a theOtherSubPaving from this set.
     *  Note that, here we require that the height of the primary root cell of
     *  this set is >= the height of the primary root cell of \a theOtherSubPaving.
     */
    Void remove_from_lower( const GridTreeSubset& theOtherSubPaving );

    /*! \brief This method changes the primary cell of this GridTreeSet.
     *  We only can increase the height of the primary cell, this is why
     *  if toPCellHeight <= this->cell().height(), then nothing is done.
     */
    Void up_to_primary_cell( const Nat toPCellHeight );

  public:
    //@{
    //! \name Constructors

    /*! \brief Create a %GridTreeSet based on zero dimensions.
     *  This constructor is needed to use the Boost Serialization library.
     */
    GridTreeSet( );

    /*! \brief The new root node can only be constructed from the existing tree node.
     *  Here, \a pRootTreeNode is not copied, we simply store its pointer inside this class.
     *  Note that, \a pRootTreeNode should correspond to the root node. \a theHeight defines
     *  the height of the primary root cell corresponding to the \a pRootTreeNode node.
     */
    GridTreeSet( const Grid& theGrid, const Nat theHeight, BinaryTreeNode * pRootTreeNode );

    /*! \brief Construct a grid tree set from a single cell.
     */
    GridTreeSet( const GridCell & theGridCell );

    /*! \brief The copy constructor that actually copies all the data,
     *  including the paving tree. I.e. the new copy of the tree is created.
     */
    GridTreeSet( const GridTreeSet & theGridTreeSet );

    /*! A simple constructor that creates the [0, 1]*...*[0, 1] cell in the
     *  \a theDimension - dimensional space. Here we assume that we have a non scaling
     *  grid with no shift of the coordinates. I.e. Grid._data._origin = {0, ..., 0}
     *  and Grid._data._lengths = {1, ..., 1}. If enable == true then the cell is enabled
     */
    explicit GridTreeSet( const Nat theDimension, const Bool enable = false );

    /*! \brief Construct an empty tree. The \a theBoundingBox is used to define the lattice
     *  block (in theGrid) that will correspond to the root of the paving tree \a pRootTreeNode.
     */
    explicit GridTreeSet( const Grid& theGrid, const Bool enable = false  );

    /*! \brief Construct an empty tree. The \a theLatticeBox is used to define the lattice
     *  block (in theGrid) that will correspond to the root of the paving tree \a pRootTreeNode.
     */
    explicit GridTreeSet( const Grid& theGrid, const ExactBox & theLatticeBox );

    /*! \brief Construct the paving based on the block's coordinates, defined by: \a theLeftLowerPoint
     *  and \a theRightUpperPoint. These are the coordinates in the lattice defined by theGrid.
     *  The primary cell, enclosing the given block of cells is computed automatically and the
     *  binary tree is rooted to that cell. The \a theEnabledCells Array defines the enabled/disabled
     *  0-level cells of the block (lexicographic order).
     */
    explicit GridTreeSet( const Grid& theGrid, const IndexArray theLeftLowerPoint,
                          const IndexArray theRightUpperPoint, const BooleanArray& theEnabledCells );

    /*! \brief Creates a new paving from the user data. \a theTree is an Array representation of the binary
     * tree structure, \a theEnabledCells tells whether a node is or is not a leaf, \a theHeight gives the
     * height of the primary cell which is assumed to correspond to the root node of \a theTree.
     */
    explicit GridTreeSet( const Grid& theGrid, Nat theHeight, const BooleanArray& theTree, const BooleanArray& theEnabledCells );

    //@}

    //@{
    //! \name Cloning/Copying/Assignment

    /*! \brief The copy assignment operator, which copies all the data
     *  including the paving tree if necessary.
     */
    GridTreeSet& operator=( const GridTreeSet & theGridTreeSet );

    /*! \brief Return a new dynamically-allocated copy of the %GridTreeSet.
     *  In this case, all the data is copied.
     */
    GridTreeSet* clone() const;

    //@}

    /*! \brief Destructor, removes all the dynamically allocated data, any */
    /* GridTreeSubset referencing this %GridTreeSet becomes invalid. */
    virtual ~GridTreeSet();



    //@{
    //! \name Geometric Operations
    /* PIETER: You may prefer to make inplace operations may return
     *a reference to *this, allowing chaining of operations.
     */

    /*! \brief Clears the set (makes empty set on same grid). */
    Void clear( );

    /*! \brief Adjoin (make inplace union with) a single cell. */
    Void adjoin( const GridCell& theCell );

    /*! \brief Remove a single cell. */
    Void remove( const GridCell& theCell );

    /*! \brief Adjoin (make inplace union with) another grid paving set. */
    Void adjoin( const GridTreeSubset& theOtherSubPaving );
    Void adjoin( const SubPavingInterface& theOtherSubPaving ) {
        this->adjoin(dynamic_cast<const GridTreeSubset&>(theOtherSubPaving)); }

    /*! \brief Restrict to (make inplace intersection with) another grid paving set. */
    Void restrict( const GridTreeSubset& theOtherSubPaving );
    Void restrict( const SubPavingInterface& theOtherSubPaving ) {
        this->restrict(dynamic_cast<const GridTreeSubset&>(theOtherSubPaving)); }

    /*! \brief Remove cells in another grid paving set. */
    Void remove( const GridTreeSubset& theOtherSubPaving );
    Void remove( const SubPavingInterface& theOtherSubPaving ) {
        this->remove(dynamic_cast<const GridTreeSubset&>(theOtherSubPaving)); }

    /*! \brief Restrict to cells rooted to the primary cell with the height (at most) \a theHeight. */
    Void restrict_to_height( const Nat theHeight );

    /*! /brief Creates an over approximation for the \a theBox on \a theGrid. \a theBox
     * is in the original space coordinates. We compute the over approximation as the
     * smallest primary cell on the Grid, such that it contains \a theBox (after it's
     * mapping on \a theGrid )
     */
    GridCell smallest_enclosing_primary_cell(const UpperBox& theBox) const;
    //@}

    //@{
    //! \name Geometric Operations

    /*! \brief Join (make union of) two grid paving sets. */
    friend GridTreeSet join( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 );

    /*! \brief The intersection of two grid paving sets. Points only lying on the
     *  intersection of the boundaries of the two sets are not included in the result.
     */
    friend GridTreeSet intersection( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 );

    /*! \brief The difference of two grid paving sets. (Results in theSet1 minus theSet2) */
    friend GridTreeSet difference( const GridTreeSubset& theSet1, const GridTreeSubset& theSet2 );

    //@}

    //@{
    //! \name Geometric Approximation

    /*! \brief Adjoin an over approximation to box, computing to the given depth:
     *   \a numSubdivInDim -- defines, how many subdivisions in each dimension from the level of
     *   the zero cell we should make to get the proper cells for outer approximating \a theSet.
     *   \pre The box must have nonempty interior.
     */
    Void adjoin_over_approximation( const ExactBox& theBox, const Nat numSubdivInDim );

    /*! \brief Adjoin an outer approximation to a given set, computing to the given depth.
     *  This method computes an outer approximation for the set \a theSet on the grid \a theGrid.
     *  Note that, the depth is the total number of subdivisions (in all dimensions) of the unit
     *  cell of the grid. This method does the followig:
     * 1. Computes the smallest Primary cell enclosing \a theSet
     * 2. Allocates the paving for this cell
     * 3. Minces the paving to the level: depth + \<the primary cell height\>
     * 4. Iterates through the enabled leaf nodes of the paving (all the nodes are initially enabled)
     * 5. Disables the cells that are disjoint with the \a theSet
     */
    Void adjoin_outer_approximation( const CompactSetInterface& theSet, const Nat numSubdivInDim );
    Void adjoin_outer_approximation( const UpperBox& theBox, const Nat numSubdivInDim );

    /*! \brief Adjoin a lower approximation to a given set, computing to the given height and depth:
     *   \a numSubdivInDim -- defines, how many subdivisions in each dimension from the level of the
     *   zero cell we should make to get the proper cells for outer approximating \a theSet.
     *   A lower approximation comprises all cells intersecting a given set.
     */
    Void adjoin_lower_approximation( const OvertSetInterface& theSet, const Nat height, const Nat numSubdivInDim );

    /*! \brief Adjoin a lower approximation to a given set restricted to the given bounding box,
     *   computing to the given depth: \a numSubdivInDim -- defines, how many subdivisions in each
     *   dimension from the level of the zero cell we should make to get the proper cells for outer
     *   approximating \a theSet. A lower approximation comprises all cells intersecting a given set.
     */
    Void adjoin_lower_approximation( const OvertSetInterface& theSet, const ExactBox& bounding_box, const Nat numSubdivInDim );

    /*! \brief Adjoin a lower approximation to a given set, computing to the given depth
     *   \a numSubdivInDim -- defines, how many subdivisions in each dimension from the level of the
     *   zero cell we should make to get the proper cells for outer approximating \a theSet.
     *   A lower approximation comprises all cells intersecting a given set.
     */
    Void adjoin_lower_approximation( const LocatedSetInterface& theSet, const Nat numSubdivInDim );

    /*! \brief Adjoin an inner approximation to a given set, computing to the given height and depth:
     *   \a numSubdivInDim -- defines, how many subdivisions in each dimension from the level of the
     *   zero cell we should make to get the proper cells for outer approximating \a theSet.
     *   An inner approximation comprises all cells that are sub-cells of the given set.
     */
    Void adjoin_inner_approximation( const OpenSetInterface& theSet, const Nat height, const Nat numSubdivInDim );

    /*! \brief Adjoin an inner approximation to a given set restricted to the given bounding box,
     *   computing to the given depth: \a numSubdivInDim -- defines, how many subdivisions in each
     *   dimension from the level of the zero cell we should make to get the proper cells for outer
     *   approximating \a theSet. An inner approximation comprises all cells that are sub-cells of
     *   the given set.
     */
    Void adjoin_inner_approximation( const OpenSetInterface& theSet, const ExactBox& bounding_box, const Nat numSubdivInDim );

    Void adjoin_inner_approximation( const LowerBox& theBox, const Nat numSubdivInDim );
    //@}

    //@{
    //! \name Input/output routines.
    //@}

};

/*! \brief This class represents a cursor/Iterator that can be used to traverse a subtree. */
class GridTreeCursor {
  private:
    //The size with which the stack size will be incremented
    static const Nat STACK_SIZE_INCREMENT = 256;

    //The index of the top stack element (within the stack Array)
    //WARNING: the initial value should be -1 to indicate that
    //there are no elements in the stack
    Int _currentStackIndex;

    /* The nodes traversed to the current location. */
    Array<BinaryTreeNode*> _theStack;

    /* The subpaving to cursor on */
    const GridTreeSubset * _pSubPaving;

    GridCell _theCurrentGridCell;

    /*! \brief Push the node into the atack */
    Void push( BinaryTreeNode* pLatestNode );

    /*! \brief Pop the node from the atack, return NULL is the stack is empty */
    BinaryTreeNode* pop( );

    /*! \brief Check is the stack contains just one element, i.e. we are at the root */
    Bool is_at_the_root() const;

    /*! \brief this method is supposed to update the _theCurrentGridCell value, when the cursos moves
     * if left_or_right == false then we go left, if left_or_right == true then right
     * if left_or_right == indeterminate then we are going one level up.
     */
    Void updateTheCurrentGridCell( Kleenean left_or_right );

    friend OutputStream& operator<<(OutputStream& os, const GridTreeCursor& theGridTreeCursor);

  public:
    /*! \brief Default constructor constructs an invalid cursor. */
    GridTreeCursor();

    /*! \brief The constructor that accepts the subpaving to cursor on */
    GridTreeCursor(const GridTreeSubset * pSubPaving);

    /*! \brief The simple copy constructor, this constructor copies the pointer to
     *  GridTreeSubset and thus the internal stack information remains valid.
     */
    GridTreeCursor(const GridTreeCursor & otherCursor);

    /*! \brief This is an assignment operator that acts similar to the copy
     *  constructor. Here we copy the pointer to underlying GridTreeSubset.
     */
    GridTreeCursor & operator=(const GridTreeCursor & otherCursor);

    /*! This destructor does not deallocate the enclosed sub paving. */
    ~GridTreeCursor();

    /*! \brief Test if the current node is enabled. */
    Bool is_enabled() const;

    /*! \brief Test if the current node is disabled*/
    Bool is_disabled() const;

    /*! \brief Test if the current node is a leaf. */
    Bool is_leaf() const;

    /*! \brief Test if the current node is the root of the subtree.
     *  Returns true if the stack has only one element. */
    Bool is_root() const;

    /*! \brief Returns true if the current node is the left child of the parent node */
    Bool is_left_child() const;

    /*! \brief Returns true if the current node is the right child of the parent node */
    Bool is_right_child() const;

    //@{
    //! \name Leaf Operations

    /*! \brief Markes the leaf node as enabled, otherwise they through \a NotALeafNodeEsception */
    Void set_enabled() const;

    /*! \brief Markes the leaf node as disabled, otherwise they through \a NotALeafNodeEsception */
    Void set_disabled() const;

    //@}

    /*! \brief Move to the parent node. Throws an NotAllowedMoveException if the current node is the root.
     *   Returns a reference to itself. */
    GridTreeCursor& move_up();

    /*! \brief Move to the left child node. Throws an NotAllowedMoveException if the current node is a leaf. */
    GridTreeCursor& move_left();

    /*! \brief Move to the right child node. Throws an NotAllowedMoveException if the current node is a leaf. */
    GridTreeCursor& move_right();

    /*! \brief Move to a child node. Throws an NotAllowedMoveException if the current node is a leaf.
     * if left_or_right == false then we go left, if left_or_right == true then right
     */
    GridTreeCursor& move(Bool left_or_right);

    /*! \brief Convert to a GridCell. */
    const GridCell& cell() const;

    /*! \brief Allows to test if the two cursors are equal, this is determined by
     * the fact that they point to the same binary-tree node.
     * NOTE: if _currentStackIndex < 0 for at least one of the cursors then the
     * result is always false.
     */
    Bool operator==(const GridTreeCursor& anotherGridTreeCursor) const;

    /*! \brief The dereferencing operator which returns a
     * reference to the GridTreeSubset for the current node.
     */
    GridTreeSubset operator*();

    /*! \brief The dereferencing operator which returns a constant
     * reference to the GridTreeSubset for the current node.
     */
    const GridTreeSubset operator*() const;
};

/*! \brief This class allows to iterate through the enabled leaf nodes of GridTreeSubset.
 * The return objects for this Iterator are constant GridCells.
 */
class GridTreeConstIterator
    : public boost::iterator_facade< GridTreeConstIterator, GridCell const, boost::forward_traversal_tag >
    , public virtual ForwardConstantIteratorInterface<GridCell>
{
  private:
    /*! \brief When set to true indicates that this is the "end Iterator" */
    Bool _is_in_end_state;

    /*! \brief the cursor object that is created based on the given sub paving*/
    GridTreeCursor _pGridTreeCursor;

    friend class boost::iterator_core_access;

    //@{
    //! \name Iterator Specific

    Void increment();

    /*! \brief Returns true if:
     * both iterators are in the "end Iterator" state
     * both iterators are NOT in the "end Iterator" state
     * and they point to the same node of the same sub paving
     */
    Bool equal( GridTreeConstIterator const & theOtherIterator) const;

    GridCell const& dereference() const;

    //@}

    //@{
    //! \name Local

    /*! \brief Allows to navigate to the first (\a firstLast==true ),
     * last (\a firstLast==false) enabled leaf of the sub paving
     * Returns true if the node was successfully found. If nothing is
     * found then the cursor should be in the root node again.
     */
    Bool navigate_to(Bool firstLast);

    /*! \brief A recursive search function, that looks for the next enabled leaf in the tree.
     *  The search is performed from left to right. (Is used for forward iteration)
     */
    Void find_next_enabled_leaf();

    //@}

  public:
    //@{
    //! \name Constructors

    /*! \brief Default constructor constructs an invalid Iterator.
     *  \internal Note that a default constructor is required for compliance with the
     *  STL Trivial Iterator concept. */
    GridTreeConstIterator();

    /*! \brief The constructor that accepts the subpacing \a pSubPaving to iterate on
     * The paramerter \a firstLastNone indicatges whether we want to position the Iterator
     * on the first enabled leaf node (firstLastNone == true) or the last one (firstLastNone == false)
     * or we are constructing the "end Iterator" that does not point anywhere.
     */
    explicit GridTreeConstIterator( const GridTreeSubset * pSubPaving, const Kleenean firstLastNone );

    /*! \brief The copy constructor that copies the current state by copying the
     *   underlying GridTreeCursor. The latter is copied by it's copy constructor.
     */
    GridTreeConstIterator( const GridTreeConstIterator& theGridPavingIter );

    //@}

    /*! \brief An assignment operator that copies the current state by copying the
     *  underlying GridTreeCursor. The latter is copied by it's assignment operator.
     */
    GridTreeConstIterator & operator=( const GridTreeConstIterator& theGridPavingIter );


    /*! This destructor only deallocates the GridTreeCursor, note that the
     * latter one does not deallocate the enclosed sub paving.
     */
    ~GridTreeConstIterator();

    //@{
    //! \name Get the cursor of the Iterator

    //This cursor is only needed to get access to enable/disable node functionality
    GridTreeCursor const& cursor() const;

    //@}
  private:
    virtual GridTreeConstIterator* clone() const { return new GridTreeConstIterator(*this); }
    virtual Bool equals( ForwardConstantIteratorInterface<GridCell> const & theOtherIterator) const {
        GridTreeConstIterator const* theOtherIteratorPointer = dynamic_cast<GridTreeConstIterator const*>(&theOtherIterator);
        return theOtherIteratorPointer && (this->equal(*theOtherIteratorPointer)); }
    virtual Void write(OutputStream& os) const { os << "GridTreeConstIterator(" << this->cursor() << ")"; }

};

/****************************************************************************************************/
/***************************************Inline functions*********************************************/
/****************************************************************************************************/

/****************************************BinaryTreeNode**********************************************/

inline Void BinaryTreeNode::init( Kleenean isEnabled, BinaryTreeNode* pLeftNode, BinaryTreeNode* pRightNode ){
    _isEnabled = isEnabled;
    _pLeftNode = pLeftNode;
    _pRightNode = pRightNode;
}

inline BinaryTreeNode::BinaryTreeNode(const Kleenean isEnabled){
    init( isEnabled, NULL, NULL );
}

inline BinaryTreeNode::BinaryTreeNode(const BinaryTreeNode& theTreeNode){
    if( (theTreeNode._pLeftNode) != NULL && ( theTreeNode._pRightNode != NULL ) ){
        init( theTreeNode._isEnabled, new BinaryTreeNode( *theTreeNode._pLeftNode ), new BinaryTreeNode( *theTreeNode._pRightNode ) );
    }else{
        //NOTE: We do not allow for nodes where one leaf is NULL and another is not
        init( theTreeNode._isEnabled, NULL, NULL );
    }
}

inline BinaryTreeNode::BinaryTreeNode( const BooleanArray& theTree, const BooleanArray& theEnabledCells ) {
    //Make default initialization
    init( false, NULL, NULL ) ;

    //If the tree is not empry and there are enabled leafs then do the thing
    if( ( theTree.size() ) > 0 && ( theEnabledCells.size() > 0 ) ) {
        Nat arr_index = 0, leaf_counter = 0;
        restore_node( this, arr_index, leaf_counter, theTree, theEnabledCells );
    } else {
        //Otherwise settle with one disabled node
        this->set_disabled();
    }
}

inline BinaryTreeNode::~BinaryTreeNode(){
    if( _pLeftNode != NULL ) {
        delete _pLeftNode;
    }
    if( _pRightNode != NULL ) {
        delete _pRightNode;
    }
}

inline Bool BinaryTreeNode::is_enabled() const{
    return definitely(_isEnabled);
}

inline Bool BinaryTreeNode::is_disabled() const{
    return ! possibly(_isEnabled) ;
}

inline Bool BinaryTreeNode::is_leaf() const{
    return (_pLeftNode == NULL) && (_pRightNode == NULL);
}

inline BinaryTreeNode * BinaryTreeNode::child_node(Bool left_or_right) const {
    return left_or_right ? _pRightNode : _pLeftNode;
}

inline BinaryTreeNode * BinaryTreeNode::left_node() const {
    return _pLeftNode;
}

inline BinaryTreeNode * BinaryTreeNode::right_node() const {
    return _pRightNode;
}

inline Void BinaryTreeNode::set(Bool enabled_or_disabled) {
    _isEnabled = enabled_or_disabled;
    if( _pLeftNode != NULL ) {
        delete _pLeftNode;
    }
    if( _pRightNode != NULL ) {
        delete _pRightNode;
    }
}

inline Void BinaryTreeNode::set_enabled() {
    if ( is_leaf() ) {
        _isEnabled = true;
    } else {
        throw NotALeafNodeException(ARIADNE_PRETTY_FUNCTION);
    }
}

inline Void BinaryTreeNode::copy_from( const BinaryTreeNode * pOtherNode ){
    if( pOtherNode != NULL ){
        _isEnabled = pOtherNode->_isEnabled;
        if( _pLeftNode != NULL){ delete _pLeftNode; _pLeftNode = NULL; }
        if( _pRightNode != NULL){ delete _pRightNode; _pRightNode = NULL; }
        if( pOtherNode->_pLeftNode != NULL ){ _pLeftNode = new BinaryTreeNode( * (pOtherNode->_pLeftNode) ); }
        if( pOtherNode->_pRightNode != NULL ){ _pRightNode = new BinaryTreeNode( * (pOtherNode->_pRightNode) ); }
    }
}

inline Void BinaryTreeNode::set_disabled() {
    if ( is_leaf() ) {
        _isEnabled = false;
    } else {
        throw NotALeafNodeException(ARIADNE_PRETTY_FUNCTION);
    }
}

inline Void BinaryTreeNode::set_unknown() {
    if ( ! is_leaf() ) {
        _isEnabled = indeterminate;
    } else {
        throw IsALeafNodeException(ARIADNE_PRETTY_FUNCTION);
    }
}

inline Void BinaryTreeNode::make_leaf(Kleenean is_enabled ){
    _isEnabled = is_enabled;
    if( _pLeftNode != NULL ) { delete _pLeftNode; _pLeftNode= NULL; }
    if( _pRightNode != NULL ) { delete _pRightNode; _pRightNode= NULL; }
}

inline Void BinaryTreeNode::split() {
    if ( is_leaf() ) {
        _pLeftNode  = new BinaryTreeNode(_isEnabled);
        _pRightNode = new BinaryTreeNode(_isEnabled);
        set_unknown();
    }
}

inline Void BinaryTreeNode::add_enabled( const BinaryWord& path ){
    add_enabled( this, path, 0 );
}

inline Void BinaryTreeNode::mince(const Nat depth) {
    mince_node(this, depth);
}

inline Void BinaryTreeNode::recombine() {
    recombine_node(this);
}

inline string BinaryTreeNode::node_to_string() const {
    stringstream tmp_stream;
    tmp_stream << "BinaryTreeNode( isLeaf = " << is_leaf() << ", isEnabled = " << is_enabled() << ", isDisabled = " << is_disabled() << " )";
    return tmp_stream.str();
}

inline BinaryTreeNode& BinaryTreeNode::operator=( const BinaryTreeNode & otherNode ) {
    //Copy the node value
    _isEnabled = otherNode._isEnabled;

    //Deallocate memory for the children, if any
    if( _pLeftNode != NULL ){ delete _pLeftNode; _pLeftNode = NULL; }
    if( _pRightNode != NULL ){ delete _pRightNode; _pRightNode = NULL; }

    //Copy the children trees from the otherNode
    if( otherNode._pLeftNode != NULL ) { _pLeftNode = new BinaryTreeNode( * ( otherNode._pLeftNode ) ); }
    if( otherNode._pRightNode != NULL ) { _pRightNode = new BinaryTreeNode( * ( otherNode._pRightNode ) ); }

    return *this;
}

/********************************************GridTreeCursor***************************************/

inline GridTreeCursor::GridTreeCursor(  ) :
    _currentStackIndex(-1), _pSubPaving(0), _theCurrentGridCell( Grid(), 0, BinaryWord() ) {
}

inline GridTreeCursor::GridTreeCursor(const GridTreeCursor & otherCursor) :
    _currentStackIndex(otherCursor._currentStackIndex), _theStack(otherCursor._theStack),
    _pSubPaving(otherCursor._pSubPaving), _theCurrentGridCell(otherCursor._theCurrentGridCell){
}

inline GridTreeCursor& GridTreeCursor::operator=(const GridTreeCursor & otherCursor) {
    _currentStackIndex = otherCursor._currentStackIndex;
    _theStack = otherCursor._theStack;
    _pSubPaving = otherCursor._pSubPaving;
    _theCurrentGridCell = otherCursor._theCurrentGridCell;

    return *this;
}

inline GridTreeCursor::GridTreeCursor(const GridTreeSubset * pSubPaving) :
    _currentStackIndex(-1), _pSubPaving(pSubPaving), _theCurrentGridCell( pSubPaving->cell() ) {
    //Remember that GridTreeSubset contains GridCell

    //Add the current node to the stack, since we are in it

    //IVAN S ZAPREEV:
    //NOTE: There is no need in allocating _theStack elements,
    //it is all done in the push method
    push(pSubPaving->_pRootTreeNode);
}

inline GridTreeCursor::~GridTreeCursor() {
    //IVAN S ZAPREEV:
    //WARNING: The subpaving should not be deallocated here!
    //There are no other data feels that we need to deallocate
    _pSubPaving = NULL;
}

inline Void GridTreeCursor::push( BinaryTreeNode* pLatestNode ){
    //If we are out of free space, increase the Array's capacity
    if( static_cast<Nat>(_currentStackIndex +1) == ( _theStack.size()  ) ){
        _theStack.resize( _theStack.size() + STACK_SIZE_INCREMENT );
    }

    //The index for the tree node to be added
    _currentStackIndex += 1;

    //Put the tree node pointer into the stack
    _theStack[_currentStackIndex] = pLatestNode;

}

inline BinaryTreeNode* GridTreeCursor::pop( ){
    BinaryTreeNode* pLastNode = NULL;

    //If there are non-root nodes in the stack
    if( _currentStackIndex > 0 ){
        //Return the stack element at _currentStackIndex,
        //then decrement the current element's index.
        return _theStack[ _currentStackIndex-- ];
    }

    return pLastNode;
}

inline Bool GridTreeCursor::is_at_the_root() const {
    //If we are pointing at the zero cell of the Array then it
    //means that we are in the root of the tree
    return (_currentStackIndex == 0);
}

inline Bool GridTreeCursor::is_enabled() const {
    return _theStack[ _currentStackIndex ]->is_enabled();
}

inline Bool GridTreeCursor::is_disabled() const {
    return _theStack[ _currentStackIndex ]->is_disabled();
}

inline Bool GridTreeCursor::is_leaf() const {
    return _theStack[ _currentStackIndex ]->is_leaf();
}

inline Bool GridTreeCursor::is_root() const {
    return  is_at_the_root();
}

inline Void GridTreeCursor::set_enabled() const {
    return _theStack[ _currentStackIndex ]->set_enabled();
}

inline Void GridTreeCursor::set_disabled() const {
    return _theStack[ _currentStackIndex ]->set_disabled();
}

inline Bool GridTreeCursor::is_left_child() const{
    //If there is a parent node and the given node is it's left child
    return ( _currentStackIndex > 0 ) && ( _theStack[ _currentStackIndex - 1 ]->left_node() == _theStack[ _currentStackIndex ] );
}

inline Bool GridTreeCursor::is_right_child() const{
    //If there is a parent node and the given node is it's right child
    return ( _currentStackIndex > 0 ) && ( _theStack[ _currentStackIndex - 1 ]->right_node() == _theStack[ _currentStackIndex ] );
}

inline GridTreeCursor& GridTreeCursor::move_up() {
    if( ! is_root() ){
        //Remove the current node from the stack
        pop();
        //Recompute the PavingGridCell
        updateTheCurrentGridCell( indeterminate );
        //Return the object back
        return ( * this);
    } else {
        throw NotAllowedMoveException(ARIADNE_PRETTY_FUNCTION);
    }
}

inline GridTreeCursor& GridTreeCursor::move_left() {
    return move(false);
}

inline GridTreeCursor& GridTreeCursor::move_right() {
    return move(true);
}

inline Void GridTreeCursor::updateTheCurrentGridCell( Kleenean left_or_right ){
    if( is_determinate(left_or_right) ){
        //Here left_or_right is either true or false, also true defines moving to the right branch
        _theCurrentGridCell._theWord.push_back( definitely( left_or_right ) );
    } else {
        //If left_or_right is indeterminate, this means that we go "up"
        //in the tree, so we remove the last bit of the path.
        _theCurrentGridCell._theWord.pop_back();
    }
    _theCurrentGridCell = GridCell( _theCurrentGridCell._theGrid, _theCurrentGridCell._theHeight, _theCurrentGridCell._theWord );
}

inline GridTreeCursor& GridTreeCursor::move(Bool left_or_right) {
    BinaryTreeNode* pNextNode;
    if( ! is_leaf() ){
        //If we are not in the leaf node then we can go down
        if( left_or_right ){ //true moves us to the right
            pNextNode = _theStack[ _currentStackIndex ]->right_node();
        } else { //false moves us to the left
            pNextNode = _theStack[ _currentStackIndex ]->left_node();
        }
        //Put the node into the stack
        push(pNextNode);
        //Recompute the PavingGridCell
        updateTheCurrentGridCell( left_or_right );
        //Return the object back
        return ( * this);
    } else {
        throw NotAllowedMoveException(ARIADNE_PRETTY_FUNCTION);
    }
}

inline const GridCell& GridTreeCursor::cell() const {
    return _theCurrentGridCell;
}

inline Bool GridTreeCursor::operator==(const GridTreeCursor& anotherGridTreeCursor) const {
    Bool areEqual = false;
    if( (this->_currentStackIndex >=0) && (anotherGridTreeCursor._currentStackIndex >=0 ) ){
        areEqual = (this->_theStack[this->_currentStackIndex] == anotherGridTreeCursor._theStack[anotherGridTreeCursor._currentStackIndex]);
    }
    return areEqual;
}

inline GridTreeSubset GridTreeCursor::operator*() {
    //IVAN S ZAPREEV:
    //NOTE: The first three parameters define the location of the _theStack[ _currentStackIndex ]
    //node with respect to the primary cell of the GridTreeSet.
    return GridTreeSubset( _theCurrentGridCell._theGrid,
                           _theCurrentGridCell._theHeight,
                           _theCurrentGridCell._theWord,
                           _theStack[ _currentStackIndex ] );
}

inline const GridTreeSubset GridTreeCursor::operator*() const {
    return (const GridTreeSubset & ) *(* this);
}

/****************************************GridTreeConstIterator************************************/

inline GridTreeConstIterator::GridTreeConstIterator( const GridTreeConstIterator& theGridPavingIter ) :
    _is_in_end_state(theGridPavingIter._is_in_end_state), _pGridTreeCursor(theGridPavingIter._pGridTreeCursor) {
    //Copy everything
}

inline GridTreeConstIterator& GridTreeConstIterator::operator=( const GridTreeConstIterator& theGridPavingIter ){
    _is_in_end_state = theGridPavingIter._is_in_end_state;
    _pGridTreeCursor = theGridPavingIter._pGridTreeCursor;

    return *this;
}

inline GridTreeConstIterator::GridTreeConstIterator( const GridTreeSubset * pSubPaving, const Kleenean firstLastNone ):
    _pGridTreeCursor(pSubPaving) {
    if( is_determinate( firstLastNone ) ){
        //If the first/last enabled node is not found, it means that there are no elements
        //to iterate on, then we switch to the "end Iterator" state
        _is_in_end_state = ! navigate_to( definitely( firstLastNone ) );
    } else {
        //In this case if we do nothing, the cursor in the Iterator will point to the root node
        //Since this can be the only node in the tre we should add a marker that indicates that
        //the Iterator is at the end state.
        _is_in_end_state = true;
    }
}

inline Void GridTreeConstIterator::increment() {
    //If we are not done iterating
    if( ! _is_in_end_state){
        //We are at some enabled-leaf node and we want to find the next one.
        //The next node is somewhere on the right in the tree.
        find_next_enabled_leaf();
    }
}

inline GridTreeConstIterator::~GridTreeConstIterator() {
    //IVAN S ZAPREEV:
    //WARNING: There are no memory allocations in this class that have to be cleaned
}

inline Bool GridTreeConstIterator::equal( GridTreeConstIterator const & theOtherIterator) const {
    //Check if both iterators are in the "end Iterator" state
    Bool result = theOtherIterator._is_in_end_state && this->_is_in_end_state;

    //If not then check if the cursors are equal (i.e. if they point to the same binary-tree node of the same (subpaving)
    if( ! result ){
        if( theOtherIterator._is_in_end_state || this->_is_in_end_state ){
            //If at one is in the end state and the other is not then the answer is FALSE
            result = false;
        } else {
            result = theOtherIterator._pGridTreeCursor == this->_pGridTreeCursor;
        }
    }

    return result;
}

inline GridCell const& GridTreeConstIterator::dereference() const {
    return _pGridTreeCursor.cell();
}

inline GridTreeCursor const& GridTreeConstIterator::cursor() const {
    return _pGridTreeCursor;
}

/********************************************GridTreeSubset******************************************/

inline GridTreeSubset::GridTreeSubset( const Grid& theGrid, const Nat theHeight,
                                       const BinaryWord& theWord, BinaryTreeNode * pRootTreeNode ) :
                                       _pRootTreeNode(pRootTreeNode), _theGridCell(theGrid, theHeight, theWord) {
}

inline GridTreeSubset::GridTreeSubset( const GridTreeSubset &otherSubset ) : _pRootTreeNode(otherSubset._pRootTreeNode),
                                                                             _theGridCell(otherSubset._theGridCell) {
}

inline GridTreeSubset::~GridTreeSubset() {
    //IVAN S ZAPREEV:
    //WARNING: This method should have no implementation what so ever
    //All the synamically allocatged data should be destroyed from the
    //corresponding Paving object
}

inline Nat GridTreeSubset::compute_number_subdiv( Float64 theWidth, const Float64 theMaxWidth) const{
    //Compute the minimum number of subdivisions N as:
    //   minimum N : ( theWidth / 2^{N} ) <= theMaxWidth
    //This is equivalent to (because all values are positive)
    //   minimum N : theWidth / theMaxWidth <= 2^{N}
    //log is a continuous-increasing dunction, thus the latter is
    //equivalent to (since 0 < theMaxWidth, theWidth < infinite )
    //   minimum N : log( theWidth / theMaxWidth ) <= N
    //In computations then we need to take
    //   N = ceil( log( theWidth / theMaxWidth ) )

    //IVAN S ZAPREEV:
    //NOTE: We take log_approx, log_approx because if we use _up and _down versions to maximize the answer,
    //then it might become very inaccurate, and in case require more subdivisions than needed.
    //NOTE: We need base 2 logarithm, we get it by the rule: log_a(b)/log_a(c) = log_c(a)
    //PIETER COLLINS:
    //NOTE: Float64 now uses approximate operators by default
    Nat result = 0;
    if ( theWidth > theMaxWidth ){
        //result = (Nat) ceil( div_approx( log_approx( div_approx( theWidth, theMaxWidth ) ) , log_approx( R(2.0) ) ) );
        result = integer_cast<Nat>(ceil( div( log( div( theWidth, theMaxWidth ) ) , log( Float64(2.0) ) ) ) );
    }
    return result;
}

inline Bool GridTreeSubset::is_empty() const {
    return BinaryTreeNode::count_enabled_leaf_nodes( this->binary_tree() ) == 0;
}

inline SizeType GridTreeSubset::size() const {
    return BinaryTreeNode::count_enabled_leaf_nodes( this->binary_tree() );
}

inline DimensionType GridTreeSubset::dimension( ) const {
    return grid().dimension();
}


inline const Grid& GridTreeSubset::grid() const {
    return this->_theGridCell.grid();
}

inline GridCell GridTreeSubset::root_cell() const {
    return _theGridCell;
}

inline Void GridTreeSubset::set_root_cell(Bool enabled_or_disabled)  {
    this->_pRootTreeNode->set(enabled_or_disabled);
}

inline UpperBox GridTreeSubset::bounding_box() const {
    if(this->is_empty()) return ExactBox(this->dimension());

    GridTreeSet::ConstIterator iter=this->begin();
    UpperBox bbox = iter->box();

    for( ; iter!=this->end(); ++iter) {
        UpperBox cell = iter->box();
        for(Nat i = 0; i < cell.dimension(); ++i) {
            if(cell[i].lower().raw() < bbox[i].lower().raw()) bbox[i].set_lower(cell[i].lower());
            if(cell[i].upper().raw() > bbox[i].upper().raw()) bbox[i].set_upper(cell[i].upper());
        }
    }

    return bbox;

}

inline Int GridTreeSubset::zero_cell_subdivisions_to_tree_subdivisions( const Nat numSubdivInDim, const Nat primaryCellHeight,
                                                                        const Nat primaryToRootCellPathLength ) const {
    //Here we take the height of the primary cell that the subpaving's root cell is rooted to
    //This height times the number of dimensions is the number of subdivisions to make in
    //the primary cell to reach the level of the zero cell. This plus the numSubdivInDim
    //times the number of dimensions gives us the total number of the primary cell we have
    //to make. But, the given paving is has the root node different from the primary cell,
    //thus we have to subtract the length of the path from the primary cell to the root cell
    //of this subpaving, to get the proper number of subdivisions to make in the binary tree
    Int theNewDepth = ( primaryCellHeight + numSubdivInDim ) * _theGridCell.grid().dimension() - primaryToRootCellPathLength;
    //If the new depth is not positive then we already have the required number
    //of subdivisions so then nothing has to be done, so we return zero!
    return (theNewDepth > 0) ? theNewDepth : 0;
}

inline Void GridTreeSubset::mince( Nat numSubdivInDim ) {
    mince_to_tree_depth( zero_cell_subdivisions_to_tree_subdivisions( numSubdivInDim, _theGridCell.height(), _theGridCell.word().size() ) );
}

inline Void GridTreeSubset::mince_to_tree_depth( const Nat theNewDepth ) {
    _pRootTreeNode->mince( theNewDepth );
}

inline Void GridTreeSubset::recombine() {
    _pRootTreeNode->recombine();
}

inline Nat GridTreeSubset::depth() const {
    return _pRootTreeNode->depth();
}

inline Bool GridTreeSubset::operator==(const GridTreeSubset& anotherGridTreeSubset) const {
    return ( this->_theGridCell == anotherGridTreeSubset._theGridCell ) &&
        ( ( * this->_pRootTreeNode ) == ( * anotherGridTreeSubset._pRootTreeNode ) );
}

inline GridTreeSubset::ConstIterator GridTreeSubset::begin() const {
    return GridTreeSubset::ConstIterator(this, true);
}

inline GridTreeSubset::ConstIterator GridTreeSubset::end() const {
    return GridTreeSubset::ConstIterator(this, indeterminate);
}

inline const BinaryTreeNode * GridTreeSubset::binary_tree() const {
    return _pRootTreeNode;
}

inline GridTreeSubset& GridTreeSubset::operator=( const GridTreeSubset &otherSubset) {
    _pRootTreeNode = otherSubset._pRootTreeNode;
    _theGridCell = otherSubset._theGridCell;

    return *this;
}

inline GridTreeSubset* GridTreeSubset::_branch(Bool left_or_right) const {
    return new GridTreeSubset(this->branch(left_or_right));
}

inline ForwardConstantIteratorInterface<GridCell>* GridTreeSubset::_begin() const {
    return new GridTreeConstIterator(this->begin());
}

inline ForwardConstantIteratorInterface<GridCell>* GridTreeSubset::_end() const {
    return new GridTreeConstIterator(this->end());
}

inline Bool GridTreeSubset::equals(const SubPavingInterface& paving) const {
    return (*this) == dynamic_cast<const GridTreeSubset&>(paving);
}

inline Bool GridTreeSubset::subset(const SubPavingInterface& paving) const {
    return Ariadne::subset(*this, dynamic_cast<const GridTreeSubset&>(paving));
}

inline Bool GridTreeSubset::intersects(const SubPavingInterface& paving) const {
    return Ariadne::intersect(*this, dynamic_cast<const GridTreeSubset&>(paving));
}

/*********************************************GridTreeSet*********************************************/

inline GridCell GridTreeSet::smallest_enclosing_primary_cell( const UpperBox& theBox ) const {
    ARIADNE_ASSERT_MSG( this->dimension() == theBox.dimension(), "Cannot find enclosing cell for ExactBox  " << theBox << " for GridTreeSet with grid " << this->grid() );

    return GridCell::smallest_enclosing_primary_cell( theBox, this->grid() );
}

inline Void GridTreeSet::adjoin( const GridCell& theCell ) {
    ARIADNE_ASSERT_MSG( this->grid() == theCell.grid(), "Cannot adjoin GridCell with grid "<<theCell.grid()<<" to GridTreeSet with grid "<<this->grid() );
    Bool has_stopped = false;
    //Align the paving and the cell
    BinaryTreeNode* pBinaryTreeNode = align_with_cell( theCell.height(), true, false, has_stopped );

    //If we are not trying to adjoin something into an enabled sub cell of the paving
    if( ! has_stopped ){
        //Add the enabled cell to the binary tree.
        pBinaryTreeNode->add_enabled( theCell.word() );
    }
}

inline Void GridTreeSet::adjoin( const GridTreeSubset& theOtherSubPaving ) {
    ARIADNE_ASSERT_MSG( this->grid() == theOtherSubPaving.cell().grid(), "Cannot adjoin GridTreeSubset with grid "<<theOtherSubPaving.cell().grid()<<" to GridTreeSet with grid "<<this->grid() );

    Bool has_stopped = false;
    //Align the paving and the cell
    BinaryTreeNode* pBinaryTreeNode = align_with_cell( theOtherSubPaving.cell().height(), true, false, has_stopped );

    //If we are not trying to adjoin something into an enabled sub cell of the paving
    if( ! has_stopped ){
        //Now, the pBinaryTreeNode of this paving corresponds to the primary cell common with theOtherSubPaving.
        //The theOtherSubPaving's root node is defined the path theOtherSubPaving.word() which starts in the
        //node corresponding to the common primary cell.
        pBinaryTreeNode->add_enabled( theOtherSubPaving.binary_tree(), theOtherSubPaving.cell().word() );
    }
}

/**************************************FRIENDS OF BinaryTreeNode***************************************/

/*! \brief Stream insertion operator, prints out two binary arrays, one is the tree structure
 *  and the other is the true/false (enabled/disabled) values for the leaf nodes
 */
inline OutputStream& operator<<(OutputStream& output_stream, const BinaryTreeNode & binary_tree ) {
    BinaryWord tree, leaves;
    binary_tree.tree_to_binary_words( tree, leaves );
    return output_stream << "BinaryTreeNode( Tree: " << tree << ", Leaves: " << leaves << ")";
}

/****************************************FRIENDS OF GridOpenCell*******************************************/

inline OutputStream& operator<<(OutputStream& os, const GridOpenCell& theGridOpenCell ) {
    //Write the grid data to the string stream
    return os << "GridOpenCell( " << theGridOpenCell.grid() <<
        ", Primary cell height: " << theGridOpenCell.height() <<
        ", Base cell path: " << theGridOpenCell.word() <<
        ", ExactBox (closure): " << theGridOpenCell.box() << " )";
}

/****************************************FRIENDS OF GridCell*******************************************/

/*! \brief Stream insertion operator for the GridCell. */
inline OutputStream& operator<<(OutputStream& os, const GridCell& gridPavingCell){
    //Write the grid data to the string stream
    return os << "GridCell( " << gridPavingCell.grid() <<
        ", Primary cell height: " << gridPavingCell.height() <<
        ", Path to the root: " << gridPavingCell.word() <<
        ", ExactBox: " << gridPavingCell.box() << " )";
}


/*************************************FRIENDS OF GridTreeCursor*****************************************/

inline OutputStream& operator<<(OutputStream& os, const GridTreeCursor& theGridTreeCursor){
    const Int curr_stack_idx = theGridTreeCursor._currentStackIndex;
    os << "GridTreeCursor( " << theGridTreeCursor._pSubPaving <<
        ", Curr. stack index: " << curr_stack_idx  <<
        ", Stack data: [ ";
    for(Int i = 0; i <= curr_stack_idx; i++){
        os << theGridTreeCursor._theStack[i]->node_to_string() << ( ( i < curr_stack_idx ) ? "" : ", ");
    }
    return os<<" ], " << theGridTreeCursor._theCurrentGridCell << " )";
}

/*************************************FRIENDS OF GridTreeSubset*****************************************/

inline OutputStream& operator<<(OutputStream& os, const GridTreeSubset& theGridTreeSubset) {
    return os << "GridTreeSubset( Primary cell: " << theGridTreeSubset.cell() << ", " << (*theGridTreeSubset.binary_tree()) <<" )";
}

/***************************************FRIENDS OF GridTreeSet******************************************/

inline OutputStream& operator<<(OutputStream& os, const GridTreeSet& theGridTreeSet) {
    const GridTreeSubset& theGridTreeSubset = theGridTreeSet;
    return os << "GridTreeSet( " << theGridTreeSubset << " )";
}

inline GridTreeSet outer_approximation( const CompactSetInterface& theSet, const Nat numSubdivInDim ) {
    Grid theGrid( theSet.dimension() );
    return outer_approximation( theSet, theGrid, numSubdivInDim );
}

inline GridTreeSet outer_approximation( const CompactSetInterface& theSet, const Grid& theGrid, const Nat numSubdivInDim ) {
    GridTreeSet result( theGrid );
    result.adjoin_outer_approximation( theSet, numSubdivInDim );
    return result;
}

inline GridTreeSet inner_approximation( const OpenSetInterface& theSet, const Grid& theGrid, const Nat height, const Nat numSubdivInDim ) {
    GridTreeSet result( theGrid );
    result.adjoin_inner_approximation( theSet, height, numSubdivInDim );
    return result;
}


template<class BS>
GridTreeSet outer_approximation(const ListSet<BS>& theSet, const Grid& theGrid, const Nat numSubdivInDim) {
    ARIADNE_ASSERT_MSG( theSet.dimension()==theGrid.dimension(),"theSet="<<theSet<<", theGrid="<<theGrid );
    GridTreeSet result( theGrid );
    for(typename ListSet<BS>::ConstIterator iter=theSet.begin(); iter!=theSet.end(); ++iter) {
        result.adjoin_outer_approximation( *iter, numSubdivInDim );
    }
    result.recombine();
    return result;
}

#ifdef ARIADNE_ENABLE_SERIALIZATION
  template<class A> Void serialize(A& archive, Ariadne::GridTreeSet& set, const unsigned int version) {
      ARIADNE_NOT_IMPLEMENTED;
  }
#endif /* ARIADNE_ENABLE_SERIALIZATION */

} // namespace Ariadne

#endif /* ARIADNE_GRID_SET_H */

