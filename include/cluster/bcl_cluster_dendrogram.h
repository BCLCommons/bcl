// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

#ifndef BCL_CLUSTER_DENDROGRAM_H_
#define BCL_CLUSTER_DENDROGRAM_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_cluster_node.h"
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically
#if defined (__GNUC__)
  #if (__GNUC__ > 3 && __GNUC_MINOR__ > 2)
    #include <backward/hash_map>
  #else
    #include <ext/hash_map>
  #endif
  using namespace __gnu_cxx;
#elif defined (_MSC_VER)
  #include <hash_map>
  using namespace std;
  using namespace stdext;
#endif

namespace bcl
{
  namespace cluster
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class KeyCompare
    //! @brief comparison type for map, set at runtime based on constructor
    //!
    //! @see @link example_cluster_dendrogram.cpp @endlink
    //! @author alexanns, woetzen
    //! @date June 6, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_PrecisionType>
    class KeyCompare
    {
    private:
      util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > m_Comparison;

    public:
      KeyCompare< t_PrecisionType>( const util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > &COMPARISON) :
        m_Comparison( COMPARISON)
      {
      }

      bool operator()( const t_PrecisionType &LHV, const t_PrecisionType &RHV) const
      {
        return m_Comparison->operator()( LHV, RHV);
      }
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Dendrogram
    //! @brief Dendrogram class is for creating the hierarchy of cluster::Nodes which is the clustering dendrogram.
    //! Upon construction the clustering is automatically performed.
    //!
    //! @see @link example_cluster_dendrogram.cpp @endlink
    //! @author alexanns, woetzen
    //! @date June 6, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class Dendrogram :
      public util::ObjectInterface
    {
    protected:

    //////////
    // data //
    //////////

      //! "m_Linkage" is the object which determines how the distance between two clusters is calculated
      util::ShPtr< LinkageInterface< t_DataType, t_PrecisionType> > m_Linkage;

      //! "m_Data" is the data which is desired to be clustered
      util::ShPtr< storage::List< t_DataType> > m_Data;

      //! "m_Node" is the Node object which will have the results of the clustering
      Node< t_DataType, t_PrecisionType> m_Node;

      //! "m_HeightCutoff" is the threshold at which the clustering is terminated
      t_PrecisionType m_HeightCutoff;

      //! ShPtr to a BinaryFunctionInterface "m_BinaryPredicate" to defines how two t_DataType are compared
      util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > m_BinaryPredicate;

      //! the threshold for combining clusters during preclustering
      t_PrecisionType m_PreClusterThreshold;

      //! used to give every node that is clustered a unique identification number
      size_t m_Identifier;

      //! similarity threshold to remove nodes at the base of the dendrogram
      t_PrecisionType m_SimilarityThreshold;

      //! the girth that is given to base nodes (i.e. nodes which are the base for building up the rest)
      static const t_PrecisionType &GetInitialNodeGirth()
      {
        return util::GetUndefined< t_PrecisionType>();
      }

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Dendrogram() :
        m_Linkage(),
        m_Data(),
        m_Node(),
        m_HeightCutoff(),
        m_BinaryPredicate(),
        m_PreClusterThreshold(),
        m_Identifier( 0),
        m_SimilarityThreshold( util::GetUndefined< t_PrecisionType>())
      {
      }

      //! @brief constructor taking ShPtr to t_DataType and a ShPtr to DistanceInterface
      //! @param DISTANCE is a ShPtr to a distance interface which will define "m_Linkage"
      //! @param DATA is the ShPtrList to t_DataType which will be clustered
      //! @param BINARY_PREDICATE defines how two t_DataType are compared
      //! @param PRECLUSTER_THRESHOLD
      Dendrogram
      (
        const util::ShPtr< LinkageInterface< t_DataType, t_PrecisionType> > &DISTANCE,
        const util::ShPtr< storage::List< t_DataType> > &DATA,
        const util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> &BINARY_PREDICATE =
          ( **math::Comparisons< t_PrecisionType>::GetEnums().e_Less),
        const t_PrecisionType PRECLUSTER_THRESHOLD = util::GetUndefined< t_PrecisionType>(),
        const t_PrecisionType SIMILARITY_THRESHOLD = util::GetUndefined< t_PrecisionType>()
      ) :
        m_Linkage( DISTANCE),
        m_Data( DATA),
        m_Node(),
        m_HeightCutoff(),
        m_BinaryPredicate( BINARY_PREDICATE.Clone()),
        m_PreClusterThreshold( PRECLUSTER_THRESHOLD),
        m_Identifier( 0),
        m_SimilarityThreshold( SIMILARITY_THRESHOLD)
      {
        util::Stopwatch clustering_timer( " dendrogram clustering ");

        m_Node = operator()();
      }

      //! @brief constructor taking ShPtr to t_DataType and a ShPtr to DistanceInterface and a t_PrecisionType
      //! @param DISTANCE is a ShPtr to a distance interface which will define "m_Linkage"
      //! @param DATA is the ShPtrList to t_DataType which will be clustered
      //! @param HEIGHT_CUTOFF is the threshold at which the clustering will be terminated
      //! @param BINARY_PREDICATE defines how two t_DataType are compared
      Dendrogram
      (
        const util::ShPtr< LinkageInterface< t_DataType, t_PrecisionType> > &DISTANCE,
        const util::ShPtr< storage::List< t_DataType> > &DATA,
        const t_PrecisionType HEIGHT_CUTOFF,
        const util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> &BINARY_PREDICATE =
          ( **math::Comparisons< t_PrecisionType>::GetEnums().e_Less),
        const t_PrecisionType PRECLUSTER_THRESHOLD = util::GetUndefined< t_PrecisionType>(),
        const t_PrecisionType SIMILARITY_THRESHOLD = util::GetUndefined< t_PrecisionType>()
      ) :
        m_Linkage( DISTANCE),
        m_Data( DATA),
        m_Node(),
        m_HeightCutoff( HEIGHT_CUTOFF),
        m_BinaryPredicate( BINARY_PREDICATE.Clone()),
        m_PreClusterThreshold( PRECLUSTER_THRESHOLD),
        m_Identifier( 0),
        m_SimilarityThreshold( SIMILARITY_THRESHOLD)
      {
        util::Stopwatch clustering_timer
        (
          " dendrogram clustering ", util::Time( std::numeric_limits< size_t>::max(), 0), util::Message::e_Critical
        );
        m_Node = operator()( HEIGHT_CUTOFF);
      }

      //! @brief Clone function
      //! @return pointer to new LinkageSingle< t_DataType>
      Dendrogram< t_DataType, t_PrecisionType> *Clone() const
      {
        return new Dendrogram< t_DataType, t_PrecisionType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief GetNode returns a const reference to "m_Node"
      //! @return returns an Node which is "m_Node"
      const Node< t_DataType, t_PrecisionType> &GetNode() const
      {
        return m_Node;
      }

      //! @brief GetHeightCutoff gives "m_HeightCutoff"
      //! @return returns the height cutoff of the dendrogram which is "m_HeightCutoff"
      t_PrecisionType GetHeightCutoff() const
      {
        return m_HeightCutoff;
      }

      //! @brief GetData gives "m_Data" which are the objects being clustered
      //! @return returns the "m_Data" which are the objects being clustered
      const util::ShPtr< storage::List< t_DataType> > &GetData() const
      {
        return m_Data;
      }

      //! @brief GetLinkage gives "m_Linkage" which determines how the distance between two clusters is calculated
      //! @return "m_Linkage"
      const util::ShPtr< LinkageInterface< t_DataType, t_PrecisionType> > &GetLinkage() const
      {
        return m_Linkage;
      }

      //! @brief GetBinaryPredicate gives "m_BinaryPredicate" which defines how two t_DataType are compared
      //! @return "m_BinaryPredicate"
      const util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > &GetBinaryPredicate() const
      {
        return m_BinaryPredicate;
      }

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

    protected:

      //! @brief operator() creates the hierarchy of Nodes
      //! @return returns Node which has the results of the clustering
      Node< t_DataType, t_PrecisionType> operator()()
      {
        //! create Node "node" with the initial base Nodes of clustering
        Node< t_DataType, t_PrecisionType> node( InitializeNode()); //< initialize with Nodes made by "m_Data"

        // true if "" is defined - indicates user wants Preclustering to occur
        if( util::IsDefined( m_PreClusterThreshold))
        {
          // precluster objects which are closely related
          PreCluster( node);
        }

        // create hash_map "linkage_involvements" this holds the addresses of a pair of nodes as two keys
        // When given two addresses (keys) it provides an iterator to the position in the map (which stores linkages
        // and the corresponding nodes) where these two nodes are
        // this allows fast lookup of a linkage given two nodes
        hash_map
        <
          size_t,
          hash_map
          <
            size_t,
            typename std::multimap
            <
              t_PrecisionType,
              std::pair
              <
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
              >,
              KeyCompare< t_PrecisionType>
            >::iterator
          >
        > linkage_involvements;

        // create KeyCompare< t_PrecisionType> "comparison" and initialize with "m_BinaryPredicate"
        const KeyCompare< t_PrecisionType> comparison( m_BinaryPredicate);

        // create multimap "linkages_nodes" and initialize with the GetInitialLinkages function
        // this stores the linkage between every node pair as a key and iterators to the nodes as values
        // This is used to quickly determine which two nodes should be combined since this map will always be sorted
        // by linkage according to KeyCompare< t_PrecisionType>, so the two nodes of the first element can be combined
        std::multimap
        <
          t_PrecisionType,
          std::pair< typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator, typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator>,
          KeyCompare< t_PrecisionType>
        > linkages_nodes( comparison);
        GetInitialLinkages( linkages_nodes, node, linkage_involvements);

        BCL_MessageStd( "initial linkages_nodes size " + util::Format()( linkages_nodes.size()));

        // go until clustering has been completed (all nodes have been combined until only two nodes are left in node)
        while( node.GetNodes().GetSize() > 2)
        {
          // set "level" to the height at which the next step in clustering occurs
          NextLevel( node, linkages_nodes, linkage_involvements);
        }

        // true if there are two nodes left in the nodes of "node"
        if( node.GetNodes().GetSize() == 2)
        {
          // calculate and set the girth of "node" which is the distance between the two nodes it contains
          node.SetGirth( m_Linkage->operator()( **node.GetNodes().Begin(), **( ++node.GetNodes().Begin())));
        }
        // if there is not two nodes left its possible possible preclustering was too ambitious or there was nothing
        // to cluster to begin with
        else
        {
          BCL_MessageCrt
          (
            "strange number of nodes left" + util::Format()( node.GetNodes().GetSize())
          );
        }

        // return "node_list" which is the list of Nodes whose girth is less than "HEIGHT_CUTOFF"
        return node;
      }

      //! @brief operator() creates the hierarchy of Nodes
      //! @param HEIGHT_CUTOFF determines the difference between two clusters that stops the agglomeration process
      //!        (inclusive for the cutoff)
      //! @return returns Node which has the results of the clustering
      Node< t_DataType, t_PrecisionType> operator()( const t_PrecisionType HEIGHT_CUTOFF)
      {
        //! create Node "node" with the initial base Nodes of clustering
        Node< t_DataType, t_PrecisionType> node( InitializeNode()); //< initialize with Nodes made by "m_Data"

        // true if "" is defined - indicates user wants Preclustering to occur
        if( util::IsDefined( m_PreClusterThreshold))
        {
          // precluster objects which are closely related
          PreCluster( node);
        }

        // create hash_map "linkage_involvements" this holds the addresses of a pair of nodes as two keys
        // When given two addresses (keys) it provides an iterator to the position in the map (which stores linkages
        // and the corresponding nodes) where these two nodes are
        // this allows fast lookup of a linkage given two nodes
        hash_map
        <
          size_t,
          hash_map
          <
            size_t,
            typename std::multimap
            <
              t_PrecisionType,
              std::pair
              <
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
              >,
              KeyCompare< t_PrecisionType>
            >::iterator
          >
        > linkage_involvements;

        // create KeyCompare< t_PrecisionType> "comparison" and initialize with "m_BinaryPredicate"
        const KeyCompare< t_PrecisionType> comparison( m_BinaryPredicate);

        // create multimap "linkages_nodes" and initialize with the GetInitialLinkages function
        // this stores the linkage between every node pair as a key and iterators to the nodes as values
        // This is used to quickly determine which two nodes should be combined since this map will always be sorted
        // by linkage according to KeyCompare< t_PrecisionType>, so the two nodes of the first element can be combined
        std::multimap
        <
          t_PrecisionType,
          std::pair< typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator, typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator>,
          KeyCompare< t_PrecisionType>
        > linkages_nodes( comparison);
        GetInitialLinkages( linkages_nodes, node, linkage_involvements);

        BCL_MessageStd( "initial linkages_nodes size " + util::Format()( linkages_nodes.size()));

        // go until the desired distance/height cutoff is reached
        while
        (
          ( m_BinaryPredicate->operator()( node.GetGirth(), HEIGHT_CUTOFF) || !util::IsDefined( node.GetGirth()))
          && node.GetNodes().GetSize() > 2
        )
        {
          // set "level" to the height at which the next step in clustering occurs
          NextLevel( node, linkages_nodes, linkage_involvements);
        }

        // true if the number of nodes left is 2 - usually will not be true since given a height cutoff the clustering
        // will usually stop before everything is clustered
        if( node.GetNodes().GetSize() == 2)
        {
          // calculate and set the girth of "node" which is the distance between the two nodes it contains
          node.SetGirth( m_Linkage->operator()( **node.GetNodes().Begin(), **( ++node.GetNodes().Begin())));
        }

        // return "node_list" which is the list of Nodes whose girth is less than "HEIGHT_CUTOFF"
        return node;
      }

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read in members
        io::Serialize::Read( m_Linkage        , ISTREAM);
        io::Serialize::Read( m_Data           , ISTREAM);
        io::Serialize::Read( m_Node           , ISTREAM);
        io::Serialize::Read( m_HeightCutoff   , ISTREAM);
        io::Serialize::Read( m_BinaryPredicate, ISTREAM);

        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_Linkage        , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Data           , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Node           , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_HeightCutoff   , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_BinaryPredicate, OSTREAM, INDENT);

        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief InitializeNodes uses "m_Data" in order to create the initial base nodes which are clustered together
      //! @return returns a ShPtrList to Nodes which is the start for the clustering
      Node< t_DataType, t_PrecisionType> InitializeNode()
      {
        util::Stopwatch clustering_timer
        (
          " initializenode ", util::Time( std::numeric_limits< size_t>::max(), 0), util::Message::e_Critical
        );
        // create ShPtrList of Nodes "node_list" which will hold the initial node of each element
        util::ShPtrList< Node< t_DataType, t_PrecisionType> > node_list;

        // iterate through "m_Data" and create a Node for each element in "m_Data" and add it to "node_list"
        for
        (
          typename storage::List< t_DataType>::iterator itr( m_Data->Begin()), itr_end( m_Data->End());
          itr != itr_end;
          ++itr
        )
        {
          // add a Node to the end of "node_list" created from the t_DataType pointed to by "itr"
          util::SiPtr< t_DataType> pointer( *itr);
          util::SiPtrList< t_DataType> member( 1, pointer);
          node_list.PushBack
          (
            // create a Node and initialize with a Node constructed from the t_DataType pointed to by "itr"
            util::ShPtr< Node< t_DataType, t_PrecisionType> >( new Node< t_DataType, t_PrecisionType>( member, GetIdentifier()))
          );
        }

        // create Node "return_node" and initialize with "node_list" and a girth of GetInitialNodeGirth()
        Node< t_DataType, t_PrecisionType> return_node( node_list, GetInitialNodeGirth(), GetIdentifier());

        // return "node" which has a Node for each t_DataType in "m_Data"
        return return_node;
      }

      //! @brief PreCluster is for going through a node and initially combining things that are closely related
      //!        This is used for quickly combining things are very closely related and would be combined anyway during
      //!        regular clustering
      //!        The difference is that this combines them in a single pass through of the nodes, while in regular
      //!        clustering it would take many iterations
      //!        This essentially uses SingleLinkage to pre-build up the clusters
      //! @param NODE this is the node where preclustering will take place
      void PreCluster( Node< t_DataType, t_PrecisionType> &NODE)
      {
        // timer to see how long preclustering takes
        util::Stopwatch precluster( "preclustering timer");

        // create map "nodes_preclustered" which will keep track for each node that is preclustered of an iterator
        // to the node into which is was engulfed and an iterator to its original individual self node
        // The iterators point to positions in the list of nodes in "NODE"
        // (engulfing nodes created during preclustering) are inserted into the list of nodes in "NODE"
        std::map
        <
          util::ShPtr< Node< t_DataType, t_PrecisionType> >,
          std::pair
          <
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator, //< containing node iterator
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator //< individual node iterator
          >
        > nodes_preclustered;

        // iterate through "NODES"
        for
        (
          typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
            itr_a( NODE.GetNodes().Begin()), itr_end( NODE.GetNodes().End());
          itr_a != itr_end;
          ++itr_a
        )
        {
          // create iterator "itr_b" and initialize with "itr_a"
          typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator itr_b( itr_a);

          // increment "itr_b" so that the second iteration through "NODES" (below) starts at one past "itr_a"
          ++itr_b;

          // second iteration through "NODES"
          for
          (
            ;
            itr_b != itr_end;
            ++itr_b
          )
          {
            // create t_PrecisionType "current_linkage" initialize with the linkage between the nodes behind "itr_a" and "itr_b"
            const t_PrecisionType current_linkage( m_Linkage->operator()( **itr_a, **itr_b));

            // true if "current_linkage" meets "m_PreClusterThreshold" i.e. the two nodes denoted by "itr_a" and "itr_b"
            // are closely related enough to combine them together into a node
            if( m_BinaryPredicate->operator()( current_linkage, m_PreClusterThreshold))
            {
              // true if neither node has been preclustered yet
              // need to combine them together into a new node that contains them both
              if
              (
                nodes_preclustered.find( *itr_a) == nodes_preclustered.end() &&
                nodes_preclustered.find( *itr_b) == nodes_preclustered.end()
              )
              {
                // create ShPtrList "combine_nodes" which will hold the two nodes
                util::ShPtrList< Node< t_DataType, t_PrecisionType> > combine_nodes;

                // add node behind "itr_a" to "combine_nodes"
                combine_nodes.PushBack( *itr_a);

                // add node behind "itr_b" to "combine_nodes"
                combine_nodes.PushBack( *itr_b);

                // create ShPtr to a Node "new_node" and initialize with "combine_nodes" and "current_linkage"
                const util::ShPtr< Node< t_DataType, t_PrecisionType> >
                  new_node( new Node< t_DataType, t_PrecisionType>( combine_nodes, current_linkage, GetIdentifier()));

                // get a reference to the list of nodes in "NODE"
                util::ShPtrList< Node< t_DataType, t_PrecisionType> > &node_nodes( NODE.GetNodes());

                // add "new_node" to the beginning of "node_nodes"
                node_nodes.PushFront
                (
                  new_node
                );

                // create const pair "containing_node_single_position_a" and initialize with the position where
                // "new_node" was inserted into "NODES" (the beginning) and the position where node_a is in the nodes
                // of NODE (i.e. "node_nodes")
                const std::pair
                <
                  typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
                  typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
                > containing_node_single_position_a( node_nodes.Begin(), itr_a);

                // add node behind "itr_a" to "nodes_preclustered" as the key and
                // "containing_node_single_position_a" as value
                nodes_preclustered[ *itr_a] = containing_node_single_position_a;

                // create const pair "containing_node_single_position_b" and initialize with the position where
                // "new_node" was inserted into "NODES" (the beginning) and the position where node_b is in the nodes
                // of NODE (i.e. "node_nodes")
                const std::pair
                <
                  typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
                  typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
                > containing_node_single_position_b( node_nodes.Begin(), itr_b);

                // add node behind "itr_b" to "nodes_preclustered" as the key and
                // "containing_node_single_position_b" as value
                nodes_preclustered[ *itr_b] = containing_node_single_position_b;
              }
              // true if node_a (i.e. the node behind "itr_a") has been preclustered before but node_b
              // (i.e. the node behind "itr_b") has not
              // just need to add node_b into the node containing node_a
              else if
              (
                nodes_preclustered.find( *itr_a) != nodes_preclustered.end() &&
                nodes_preclustered.find( *itr_b) == nodes_preclustered.end()
              )
              {
                // create reference to ShPtr "node_a" and initialize with the node that contains node behind "itr_a"
                util::ShPtr< Node< t_DataType, t_PrecisionType> > &node_a( *nodes_preclustered[ *itr_a].first);

                // set the girth of "node_a" now that it is going to also contain node_b
                node_a->SetGirth( m_Linkage->operator()( *node_a, **itr_b));

                // add the node behind "itr_b" into "node_a"
                node_a->GetNodes().PushBack( *itr_b);

                // add the members of the node behind "itr_b" into the members of "node_a"
                node_a->GetMembers().Append( ( *itr_b)->GetMembers());

                // create const pair "containing_node_single_position_b" and initialize with the position where
                // "new_node" was inserted into "NODES" (the beginning) and the position where node_b is in the nodes
                // of NODE (i.e. "node_nodes")
                const std::pair
                <
                  typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
                  typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
                > containing_node_single_position_b( nodes_preclustered[ *itr_a].first, itr_b);

                // add node behind "itr_b" to "nodes_preclustered" as the key and
                // "containing_node_single_position_b" as value
                nodes_preclustered[ *itr_b] = containing_node_single_position_b;
              }
              // true if node_a (i.e. the node behind "itr_a") has not been preclustered before but node_b
              // (i.e. the node behind "itr_b") has been
              // just need to add node_a into the node containing node_b
              else if
              (
                nodes_preclustered.find( *itr_a) == nodes_preclustered.end() &&
                nodes_preclustered.find( *itr_b) != nodes_preclustered.end()
              )
              {
                // create reference to ShPtr "node_b" and initialize with the node that contains node behind "itr_b"
                util::ShPtr< Node< t_DataType, t_PrecisionType> > &node_b( *nodes_preclustered[ *itr_b].first);

                // set the girth of "node_b" now that it is going to also contain node_a
                node_b->SetGirth( m_Linkage->operator()( *node_b, **itr_a));

                // add the node behind "itr_a" into "node_b"
                node_b->GetNodes().PushBack( *itr_a);

                // add the members of the node behind "itr_a" into the members of "node_b"
                node_b->GetMembers().Append( ( *itr_a)->GetMembers());

                // create const pair "containing_node_single_position_a" and initialize with the position where
                // "new_node" was inserted into "NODES" (the beginning) and the position where node_a is in the nodes
                // of NODE (i.e. "node_nodes")
                const std::pair
                <
                  typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
                  typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
                > containing_node_single_position_a( nodes_preclustered[ *itr_b].first, itr_a);

                // add node behind "itr_a" to "nodes_preclustered" as the key and
                // "containing_node_single_position_a" as value
                nodes_preclustered[ *itr_a] = containing_node_single_position_a;
              }
              // true if both node_a (i.e. the node behind "itr_a") and node_b (i.e. the node behind "itr_b") have both
              // been preclustered previously
              // Will move all the nodes involved with node_b (i.e.all the nodes that are in the node containing node_b)
              // into the node containing node_a
              else if
              (
                nodes_preclustered.find( *itr_a) != nodes_preclustered.end() &&
                nodes_preclustered.find( *itr_b) != nodes_preclustered.end()
              )
              {
                // create reference to ShPtr "node_b" and initialize with the node that contains node behind "itr_b"
                util::ShPtr< Node< t_DataType, t_PrecisionType> > &node_b_containing_node( *nodes_preclustered[ *itr_b].first);

                // create reference to ShPtr "node_a" and initialize with the node that contains node behind "itr_a"
                util::ShPtr< Node< t_DataType, t_PrecisionType> > &node_a_containing_node( *nodes_preclustered[ *itr_a].first);

                // true if node_a and node_b are actually already contained in the same node
                // can just skip everything and continue on to next node pair
                if( node_a_containing_node == node_b_containing_node)
                {
                  continue;
                }

                // debug message printing out the nodes to be combined
                BCL_MessageDbg
                (
                  "need to combine\n" + util::Format()( node_b_containing_node) + "\nand\n" +
                  util::Format()( node_a_containing_node)
                );

                // calculate the linkage between "node_a_containing_node" and "node_b_containing_node" and set the
                // linkage of "node_a_containing_node" to this calculated girth
                node_a_containing_node->SetGirth
                (
                  m_Linkage->operator()( *node_a_containing_node, *node_b_containing_node)
                );

                // add the nodes in "node_b_containing_node" into the nodes of of "node_a_containing_node"
                node_a_containing_node->GetNodes().Append( node_b_containing_node->GetNodes());

                // add the members of "node_b_containing_node" into the members of "node_a_containing_node"
                node_a_containing_node->GetMembers().Append( node_b_containing_node->GetMembers());

                // debug message printing out "node_a_containing_node"
                BCL_MessageDbg
                (
                  "new combined node a\n" + util::Format()( node_a_containing_node)
                );

                // need to make a copy of the containing node iterator because it is going to be reset in the for loop
                // below and need this to remove the containing node b from NODE since it is no longer needed
                // cannot remove the containing node b before the for loop because we need to iterate through it to
                // reset all the contained nodes iterators to point to containing node a
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator node_b_containing_node_itr
                (
                  nodes_preclustered[ *itr_b].first
                );

                // set all the nodes contained in containing node B to know that they are now contained in containing
                // node A
                for
                (
                  typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
                    itr( node_b_containing_node->GetNodes().Begin()), itr_end( node_b_containing_node->GetNodes().End());
                  itr != itr_end;
                  ++itr
                )
                {
                  // debug messages
                  BCL_MessageDbg( "resetting containing iterator for\n" + util::Format()( *itr));
                  BCL_MessageDbg
                  (
                    "setting containing iterator to\n" + util::Format()( *nodes_preclustered[ *itr_a].first)
                  );

                  // for the node denoted by "itr" set its containing node iterator to "itr_a"
                  nodes_preclustered[ *itr].first = nodes_preclustered[ *itr_a].first;
                }

                // remove containing node b from NODES since its members are now contained in node A containing node
                NODE.GetNodes().Remove( node_b_containing_node_itr);

                // debug message printing out the updated "node_a_containing_node"
                BCL_MessageDbg
                (
                  "new combined node\n" + util::Format()( node_a_containing_node)
                );

              }
            }
          } //< end inner iteration through NODES
        } //< end outer iteration through NODES

        // create reference to ShPtrList "node_nodes" and initialize with the nodes of "NODE"
        util::ShPtrList< Node< t_DataType, t_PrecisionType> > &node_nodes( NODE.GetNodes());

        BCL_MessageDbg
        (
          "before pruning NODE.getnodes\nnode nodes\n" + util::Format()( node_nodes)
        );

        // delete the single nodes from NODE which were preclustered - they are now contained in a containing node
        // this is done by iterating through "nodes_preclustered" and taking advantage that it stores iterators
        // to the position of the original individual single node in "node_nodes"
        for
        (
          typename std::map
          <
            util::ShPtr< Node< t_DataType, t_PrecisionType> >,
            std::pair< typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator, typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator>
          >::const_iterator itr( nodes_preclustered.begin()), itr_end( nodes_preclustered.end());
          itr != itr_end;
          ++itr
        )
        {
          // remove the single individual node from "node_nodes" using the stored iterator in "nodes_preclustered"
          // which denotes its position in "node_nodes"
          node_nodes.Remove( itr->second.second);
        }

        // messages
        BCL_MessageDbg( "preclustered node is\n" + util::Format()( NODE));
        BCL_MessageCrt
        (
          "number of preclustered node is\n" + util::Format()( nodes_preclustered.size())
        );
      }

      //! @brief GetInitialLinkages calculates all the initial linkages between all the individual nodes
      //! @param NODE is the node which contains all of the initial nodes whose pairwise linkages will be calculated
      //! @param LINKAGES_NODES multimap storing all pairwise linkages sorted by linkage
      //! @param LINKAGE_INVOLVEMENTS is a hash_map of hash_map where each key is the address of one of the nodes
      //!        in NODE and the value is an iterator to the position in the linkage map
      //!        (being created by this function) where the linkage between the two nodes was inserted
      //! @return multimap which has as keys the linkage between two nodes and as value iterators to the two nodes
      void GetInitialLinkages
      (
        std::multimap
        <
          t_PrecisionType,
          std::pair
          <
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
          >,
          KeyCompare< t_PrecisionType>
        > &LINKAGES_NODES,
        Node< t_DataType, t_PrecisionType> &NODE,
        hash_map
        <
          size_t, //< node address
          hash_map
          <
            size_t, //< node address
            // iterator to position in linkage map where the linkage between the two nodes (and iterators to
            // the two nodes) is
            typename std::multimap
            <
              t_PrecisionType,
              std::pair
              <
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
              >,
              KeyCompare< t_PrecisionType>
            >::iterator
          >
        > &LINKAGE_INVOLVEMENTS
      ) const
      {
        util::Stopwatch clustering_timer
        (
          " get_initial_linkages ", util::Time( std::numeric_limits< size_t>::max(), 0), util::Message::e_Critical
        );
        // iterate through "NODES"
        for
        (
          typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
            itr_a( NODE.GetNodes().Begin()), itr_end( NODE.GetNodes().End());
          itr_a != itr_end;
          ++itr_a
        )
        {
          // create iterator "itr_b" and initialize with "itr_a"
          typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator itr_b( itr_a);

          // increment "itr_b" so that the second iteration through "NODES" (below) starts at one past itr_a
          ++itr_b;

          // create const reference "shpr_node_a" and initialize with the node behind "itr_a"
          const util::ShPtr< Node< t_DataType, t_PrecisionType> > &shpr_node_a( *itr_a);

          // create size_t "key_node_a" and initialize with the address of "shpr_node_a"
          // this will be used as the key for node_a in "LINKAGE_INVOLVEMENTS"
          const size_t key_node_a( size_t( shpr_node_a.GetPointer()));

          // create reference to hash_map "map_a" and initialize with the hash_map that is created when "key_node_a"
          // is inserted into "LINKAGE_INVOLVEMENTS"
          hash_map
          <
            size_t,
            typename std::multimap
            <
              t_PrecisionType,
              std::pair
              <
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
              >,
              KeyCompare< t_PrecisionType>
            >::iterator
          > &map_a( LINKAGE_INVOLVEMENTS[ key_node_a]);

          // second iteration through "NODES"
          for
          (
            ;
            itr_b != itr_end;
            ++itr_b
          )
          {
            // create const std::pair "node_pair" and initialize with "itr_a" and "itr_b" which are iterators to the
            // position where node_a and node_b are in the master node list
            const std::pair
            <
              typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
              typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
            > node_pair( itr_a, itr_b);

            // create const reference to a ShPtr "shpr_node_b" and initialize with the ShPtr behind "itr_b"
            const util::ShPtr< Node< t_DataType, t_PrecisionType> > &shpr_node_b( *itr_b);

            // create const size_t "key_node_b" and initialize with the address "shpr_node_b" points to
            // this will be used as the key for node_b in "LINKAGE_INVOLVEMENTS"
            const size_t key_node_b( size_t( shpr_node_b.GetPointer()));

            // create t_PrecisionType "current_linkage" initialize with linkage between the nodes denoted by "itr_a" and "itr_b"
            const t_PrecisionType current_linkage( m_Linkage->operator()( *shpr_node_a, *shpr_node_b));

            // create multimap iterator "insert_pos" and initialize with the position where "current_linkage" and
            // "node_pair" are inserted into "intial_linkages"
            // this iterator is used in "LINKAGE_INVOLVEMENTS" so that given two nodes, their linkage and iterators
            // to the actual nodes can be quickly gotten
            const typename std::multimap
            <
              t_PrecisionType,
              std::pair
              <
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
              >,
              KeyCompare< t_PrecisionType>
            >::iterator insert_pos
            (
              // insert "current_linkage" and "node_pair" into "initial_linkages"
              LINKAGES_NODES.insert
              (
                std::pair
                <
                  t_PrecisionType,
                  std::pair
                  <
                    typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
                    typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
                  >
                >( current_linkage, node_pair)
              )
            );

            // insert address of node_b as key to "map_a" with "insert_pos" as value
            map_a[ key_node_b] = insert_pos;

            // insert "key_node_b" and "key_node_a" node pair into "LINKAGE_INVOLVEMENTS" with "insert_pos"
            // this is done so that either node_a or node_b can be passed first to "LINKAGE_INVOLVEMENTS" and the
            // value will be able to be found
            LINKAGE_INVOLVEMENTS[ key_node_b][ key_node_a] = insert_pos;
          }
        }
      }

      //! @brief NextLevel generates the next step in the clustering hierarchy
      //! @param NODE is the list of Nodes which is being clustered together
      //! @param LINKAGES_NODES multimap storing all pairwise linkages sorted by linkage
      //! @param NODE_LINKAGE_INVOLVEMENT for all node pairs stores an iterator to the corresponding position in "LINKAGES_NODES"
      //!        Allows fast lookup of a linkage given two nodes
      void NextLevel
      (
        Node< t_DataType, t_PrecisionType> &NODE,
        std::multimap
        <
          t_PrecisionType,
          std::pair
          <
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
          >,
          KeyCompare< t_PrecisionType>
        > &LINKAGES_NODES,
        hash_map
        <
          size_t,
          hash_map
          <
            size_t,
            typename std::multimap
            <
              t_PrecisionType,
              std::pair
              <
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
              >,
              KeyCompare< t_PrecisionType>
            >::iterator
          >
        > &NODE_LINKAGE_INVOLVEMENT
      )
      {
        // create iterator "linkages_nodes_begin" initialize to the beginning of "LINKAGES_NODES"
        // this is the linkage and two nodes that should be combined together since "LINKAGES_NODES" is sorted by
        // KeyCompare< t_PrecisionType>
        const typename std::multimap
        <
          t_PrecisionType,
          std::pair< typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator>,
          KeyCompare< t_PrecisionType>
        >::iterator linkages_nodes_begin( LINKAGES_NODES.begin());

        // create const reference to std::pair "minimium_pair" and initialize with the pair holding the linkage
        // and the pair of iterators to nodes that should be combined
        const std::pair
        <
          t_PrecisionType,
          std::pair
          <
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
          >
        > &minimum_pair( *linkages_nodes_begin);

        // create const t_PrecisionType "minimum_linkage" and initialize with the linkage of the two nodes that will be combined
        const t_PrecisionType minimum_linkage( minimum_pair.first);

        // create ShPtrList "combine_nodes" this will hold the two nodes to be combined and be used to create the new
        // combined node containing them both
        util::ShPtrList< Node< t_DataType, t_PrecisionType> > combine_nodes;

        // get const reference "node_pair" to the pair of iterators denoting the two nodes to be combined
        const std::pair
        <
          typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
          typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
        > &node_pair( minimum_pair.second);

        BCL_MessageDbg
        (
          "combining node with members\n" + util::Format()( ( *node_pair.first)->GetMembers())
          + "\n with node with members\n" + util::Format()( ( *node_pair.second)->GetMembers())
          + "\nthey have a linkage of " + util::Format()( minimum_linkage)
        )

        // create const reference to ShPtrList iterator "node_pair_itr_a" and initialize with the first iterator in
        // "node_pair"
        const typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator &node_pair_itr_a( node_pair.first);

        // create const reference to ShPtr to node "node_a" and initialize with the ShPtr behind "node_pair_itr_a"
        const util::ShPtr< Node< t_DataType, t_PrecisionType> > &node_a( *node_pair_itr_a);

        // add "node_a" to "combine_nodes"
        combine_nodes.PushBack( node_a);

        // create const reference to ShPtrList iterator "node_pair_itr_b" and initialize with the second iterator in
        // "node_pair"
        const typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator &node_pair_itr_b( node_pair.second);

        // create const reference to ShPtr to node "node_b" and initialize with the ShPtr behind "node_pair_itr_b"
        const util::ShPtr< Node< t_DataType, t_PrecisionType> > &node_b( *node_pair_itr_b);

        // add "node_a" to "combine_nodes"
        combine_nodes.PushBack( node_b);

        // create const ShPtr to node "new_node" and initialize with node created with "combine_nodes" and
        // "minimum_linkage"
        util::ShPtr< Node< t_DataType, t_PrecisionType> > new_node( new Node< t_DataType, t_PrecisionType>( combine_nodes, minimum_linkage, GetIdentifier()));

        // true if the similarity threshold for nodes is defined
        if( util::IsDefined( m_SimilarityThreshold))
        {
          // remove any internally similar nodes from new node
          new_node->RemoveNodesWithHighSimilarity( m_SimilarityThreshold, *m_BinaryPredicate);
        }

        // create reference to ShPtrList "node_nodes" and initialize with the nodes of "NODE"
        util::ShPtrList< Node< t_DataType, t_PrecisionType> > &node_nodes( NODE.GetNodes());

        // add "new_node" to "node_nodes"
        node_nodes.PushBack
        (
          new_node
        );

        // true if there are more than four nodes in "node_nodes" - need to update linkages
        // if there not more than four nodes i.e. if there are four, it means the clustering is finished because
        // down below "node_pair_itr_a" and "node_pair_itr_b" will be removed, leaving only two nodes in "NODE"
        if( node_nodes.GetSize() > 4)
        {
          UpdateLinkages( NODE, LINKAGES_NODES, NODE_LINKAGE_INVOLVEMENT, new_node, node_a, node_b, linkages_nodes_begin);
        }

        // remove iterator which previously pointed to the beginning of "LINKAGES_NODES" since these nodes have now
        // been combined
        LINKAGES_NODES.erase
        (
            linkages_nodes_begin
        );

        // remove the Nodes from "node_nodes" denoted by "itr_a" and "itr_b"
        node_nodes.Remove( node_pair_itr_a);
        node_nodes.Remove( node_pair_itr_b);

        // set the girth of "NODE"
        NODE.SetGirth( minimum_linkage);
      }

    private:

      //! @brief UpdateLinkages is for updating information about pairwise node linkages
      //!        It is important to do after two nodes have been combined
      //! @param NODE the master node in which all the clustering is taking place
      //! @param LINKAGES_NODES multimap storing all pairwise linkages sorted by linkage
      //! @param NODE_LINKAGE_INVOLVEMENT for all node pairs stores an iterator to the corresponding position in "LINKAGES_NODES"
      //!        Allows fast lookup of a linkage given two nodes
      //! @param NEW_NODE the new node that was just created by combining two nodes
      //! @param NEW_NODE_COMPONENT_A the first of two nodes that was used to create NEW_NODE
      //! @param NEW_NODE_COMPONENT_B the second of two nodes that was used to create NEW_NODE
      //! @param COMBINED_NODES_ITR iterator to the position in LINKAGES_NODES with the linkage and nodes for combining
      void UpdateLinkages
      (
        Node< t_DataType, t_PrecisionType> &NODE,
        std::multimap
        <
          t_PrecisionType,
          std::pair< typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator, typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator>,
          KeyCompare< t_PrecisionType>
        > &LINKAGES_NODES,
        hash_map
        <
          size_t,
          hash_map
          <
            size_t,
            typename std::multimap
            <
              t_PrecisionType,
              std::pair
              <
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
              >,
              KeyCompare< t_PrecisionType>
            >::iterator
          >
        > &NODE_LINKAGE_INVOLVEMENT,
        const util::ShPtr< Node< t_DataType, t_PrecisionType> > &NEW_NODE,
        const util::ShPtr< Node< t_DataType, t_PrecisionType> > &NEW_NODE_COMPONENT_A,
        const util::ShPtr< Node< t_DataType, t_PrecisionType> > &NEW_NODE_COMPONENT_B,
        const typename std::multimap
        <
          t_PrecisionType,
          std::pair
          <
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
          >,
          KeyCompare< t_PrecisionType>
        >::iterator &COMBINED_NODES_ITR
      )
      {
        BCL_MessageDbg( "updating linkages");

        // create size_t "new_node_component_a_key" and initialize with the address of "NEW_NODE_COMPONENT_A"
        // this will be used to access information from "NODE_LINKAGE_INVOLVEMENT"
        const size_t new_node_component_a_key( size_t( NEW_NODE_COMPONENT_A.GetPointer()));

        // create size_t "new_node_component_b_key" and initialize with the address of "NEW_NODE_COMPONENT_B"
        // this will be used to access information from "NODE_LINKAGE_INVOLVEMENT"
        const size_t new_node_component_b_key( size_t( NEW_NODE_COMPONENT_B.GetPointer()));

        // create reference to hashmap "node_a_linkage_involvements" and initialize with the hash_map that maps the
        // linkages node_a is involved in
        hash_map
        <
          size_t,
          typename std::multimap
          <
            t_PrecisionType,
            std::pair
            <
              typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
              typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
            >,
            KeyCompare< t_PrecisionType>
          >::iterator
        > &node_a_linkage_involvements( NODE_LINKAGE_INVOLVEMENT[ new_node_component_a_key]);

        // iterate through "node_a_linkage_involvements" in order to recalculate linkages between the new node
        // and nodes that previously had linkages calculated against node_a
        for
        (
          typename hash_map
          <
            size_t,
            typename std::multimap
            <
              t_PrecisionType,
              std::pair
              <
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
              >,
              KeyCompare< t_PrecisionType>
            >::iterator
          >::iterator itr_a( node_a_linkage_involvements.begin()),
            itr_a_end( node_a_linkage_involvements.end());
          itr_a != itr_a_end; ++itr_a
        )
        {
          // create const reference to multimap iterator "itr_a_multimap_itr" and initialize with the iterator
          // that is denoted by "itr_a"
          // this iterator points to one of the positions in the linkage map
          const typename std::multimap
          <
            t_PrecisionType,
            std::pair
            <
              typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
              typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
            >,
            KeyCompare< t_PrecisionType>
          >::iterator &itr_a_multimap_itr( itr_a->second);

          // true if "itr_a_multimap_itr" is "COMBINED_NODES_ITR" which is the beginning of the linkage map
          // the beginning of the linkage map always has the linkage and nodes that are being combined
          // therefore if this is true we do not need to update the linkage involvement since the two nodes denoted by
          // "itr_a_multimap_itr" are the ones that were combined
          if( COMBINED_NODES_ITR == itr_a_multimap_itr)
          {
            continue;
          }

          // create hash_map iterator "itr_b" initialize with an iterator pointing to the place (in the hash_map of
          // linkages that node_b is involved in) where the node whose linkage is currently being updated exists
          // Since every node has its linkage calculated with every other node, the node-to-be-updated should definitely
          // exist in the hash_map of linkages of node_b
          // This "itr_b" is used below just to remove the node_b and node-to-be-updated pair from "LINKAGES_NODES"
          typename hash_map
          <
            size_t,
            typename std::multimap
            <
              t_PrecisionType,
              std::pair
              <
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
              >,
              KeyCompare< t_PrecisionType>
            >::iterator
          >::iterator itr_b( NODE_LINKAGE_INVOLVEMENT[ new_node_component_b_key].find( itr_a->first));

          // make sure the node-to-be-updated could be found in the hash_map of linkages of node_b
          BCL_Assert
          (
            itr_b != NODE_LINKAGE_INVOLVEMENT[ new_node_component_b_key].end(),
            util::Format()( itr_a->first) + "\nnot found in component b involvements"
          );

          // now find out which iterator in the pair denoted by "itr_a" is not pointing to "node_a" since the node that
          // is not node a is the one that needs to have its linkage with the new node calculated.

          // create const reference to std::pair "itr_a_multimap_itr_pair" and initialize with the pair behind
          // "itr_a_multimap_itr"
          // this is the pair of node_a and the node-to-be-updated, but which is which is undetermined
          const std::pair
          <
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
          > &itr_a_multimap_itr_pair( itr_a_multimap_itr->second);

          // create iterator to ShPtrList of nodes "not_node_a_itr" and guess initialize with the first iterator in
          // "itr_a_multimap_itr_pair"
          typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator not_node_a_itr( itr_a_multimap_itr_pair.first);

          // true if the node behind "not_node_a_itr" is "NEW_NODE_COMPONENT_A" (a.k.a node_a)
          if( *not_node_a_itr == NEW_NODE_COMPONENT_A)
          {
            // need to set "not_node_a_itr" to the second iterator in "itr_a_multimap_itr_pair"
            not_node_a_itr = itr_a_multimap_itr_pair.second;
          }

          // create const reference to ShPtr to node "not_node_a_node" initialize with the node behind "not_node_a_itr"
          const util::ShPtr< Node< t_DataType, t_PrecisionType> > &not_node_a_node( *not_node_a_itr);

          // create const t_PrecisionType "current_linkage" initialize with the linkage between "NEW_NODE" and "not_node_a_node"
          const t_PrecisionType current_linkage( m_Linkage->operator()( *NEW_NODE, *not_node_a_node));

          // create std::pair "new_node_pair" and initialize with iterators to the last node in the nodes of "NODE"
          // and the "not_node_a_itr"
          // the last node in the nodes of "NODE" is the new node, so this pair contains the pair of the new node and
          // "not_node_a_itr"
          std::pair
          <
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
          > new_node_pair( --NODE.GetNodes().End(), not_node_a_itr);

          // create multimap iterator "insert_pos" and initialize with the position where "current_linkage" and its
          // corresponding value "new_node_pair" were inserted
          const typename std::multimap
          <
            t_PrecisionType,
            std::pair
            <
              typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
              typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
            >,
            KeyCompare< t_PrecisionType>
          >::iterator insert_pos
          (
            // insert "current_linkage" and "new_node_pair" as a pair into "LINKAGES_NODES"
            // so now the linkage map has the updated linkage between the new node and "not_node_a_itr"
            LINKAGES_NODES.insert
            (
              std::pair
              <
                t_PrecisionType,
                std::pair
                <
                  typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
                  typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
                >
              >( current_linkage, new_node_pair)
            )
          );

          // erase the linkage node pair element of "LINKAGES_NODES" that NEW_NODE_COMPONENT_A was involved in
          LINKAGES_NODES.erase( itr_a_multimap_itr);

          // erase the linkage node pair element of "LINKAGES_NODES" that NEW_NODE_COMPONENT_B was involved in
          LINKAGES_NODES.erase( itr_b->second);

          // create const size_t "not_node_a_node_key" and initialize with the address pointed to by "not_node_a_node"
          // this will be used as a key for one of the linkages "NEW_NODE" is involved in
          // it will also be used to add the "NEW_NODE" linkage to the list of linkages "not_node_a_node" is involved in
          const size_t not_node_a_node_key( size_t( not_node_a_node.GetPointer()));

          // create reference to hash_map "current_involvement" and initialize with the linkages that not_node_a
          // is involved in
          // not_node_a's involvement with "NEW_NODE" needs to be added to this hash_map
          hash_map
          <
            size_t,
            typename std::multimap
            <
              t_PrecisionType,
              std::pair
              <
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator,
                typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
              >,
              KeyCompare< t_PrecisionType>
            >::iterator
          > &current_involvement( NODE_LINKAGE_INVOLVEMENT[ not_node_a_node_key]);

          // create size_t "new_node_key" and initialize with the address pointed to by "NEW_NODE"
          const size_t new_node_key( size_t( NEW_NODE.GetPointer()));

          // set the value of the combination of "new_node_key" and "not_node_a_node_key" in "NODE_LINKAGE_INVOLVEMENT"
          // to be "insert_pos" so that given "NEW_NODE" and not_node_a, the linkage between them can be gotten
          NODE_LINKAGE_INVOLVEMENT[ new_node_key][ not_node_a_node_key] = insert_pos;

          // set the value of "new_node_key" in "current_involvement" to be "insert_pos" so that given not_node_a
          // and "NEW_NODE" (the opposite order as above), the linkage between them can be gotten
          current_involvement[ new_node_key] = insert_pos;

          // need to remove node_a involvement from the list of involvements of not_node_a
          current_involvement.erase( new_node_component_a_key);

          // need to remove node_b involvement from the list of linkage involvements of not_node_a
          current_involvement.erase( new_node_component_b_key);

        }

        // need to erase node_a from "NODE_LINKAGE_INVOLVEMENT" since it no longer exists except as part of NEW_NODE
        NODE_LINKAGE_INVOLVEMENT.erase( new_node_component_a_key);

        // need to erase node_b from "NODE_LINKAGE_INVOLVEMENT" since it no longer exists except as part of NEW_NODE
        NODE_LINKAGE_INVOLVEMENT.erase( new_node_component_b_key);
      } //< end UpdateLinkages

      //! @brief gives a unique identifying number for each node in the dendrogram
      //! @return a size_t which has not be returned before
      size_t GetIdentifier()
      {
        // increase the variable and return the increased value
        return ++m_Identifier;
      }

    }; // Dendrogram

    // instantiate s_Instance
    template< typename t_DataType, typename t_PrecisionType>
    const util::SiPtr< const util::ObjectInterface> Dendrogram< t_DataType, t_PrecisionType>::s_Instance
    (
      GetObjectInstances().AddInstance( new Dendrogram< t_DataType, t_PrecisionType>())
    );

    //! @brief operator == checks if two Dendrograms are the same
    //! @param DENDROGRAM_A first Dendrogram
    //! @param DENDROGRAM_B second Dendrogram
    //! @return true, if Dendrograms are identical
    template< typename t_DataType, typename t_PrecisionType>
    inline bool operator ==
    (
      const Dendrogram< t_DataType, t_PrecisionType> &DENDROGRAM_A,
      const Dendrogram< t_DataType, t_PrecisionType> &DENDROGRAM_B
    )
    {
      return DENDROGRAM_A.GetLinkage() == DENDROGRAM_B.GetLinkage()                 &&
             DENDROGRAM_A.GetNode() == DENDROGRAM_B.GetNode()                       &&
             DENDROGRAM_A.GetHeightCutoff() == DENDROGRAM_B.GetHeightCutoff()       &&
             DENDROGRAM_A.GetData() == DENDROGRAM_B.GetData()                       &&
             DENDROGRAM_A.GetBinaryPredicate() == DENDROGRAM_B.GetBinaryPredicate();
    }

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_DENDROGRAM_H_
