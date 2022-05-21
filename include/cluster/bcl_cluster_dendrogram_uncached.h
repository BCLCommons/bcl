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

#ifndef BCL_CLUSTER_DENDROGRAM_UNCACHED_H_
#define BCL_CLUSTER_DENDROGRAM_UNCACHED_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_cluster_dendrogram.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DendrogramUncached
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_cluster_dendrogram_uncached.cpp @endlink
    //! @author alexanns
    //! @date Jun 2, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class DendrogramUncached :
      public Dendrogram< t_DataType, t_PrecisionType>
    {

    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DendrogramUncached< t_DataType, t_PrecisionType>() :
        Dendrogram< t_DataType, t_PrecisionType>()
      {
      }

      //! @brief constructor taking ShPtr to t_DataType and a ShPtr to DistanceInterface
      //! @param DISTANCE is a ShPtr to a distance interface which will define "m_Linkage"
      //! @param DATA is the ShPtrList to t_DataType which will be clustered
      //! @param BINARY_PREDICATE defines how two t_DataType are compared
      //! @param PRECLUSTER_THRESHOLD
      DendrogramUncached< t_DataType, t_PrecisionType>
      (
        const util::ShPtr< LinkageInterface< t_DataType, t_PrecisionType> > &DISTANCE,
        const util::ShPtr< storage::List< t_DataType> > &DATA,
        const util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> &BINARY_PREDICATE =
          ( **math::Comparisons< t_PrecisionType>::GetEnums().e_Less),
        const t_PrecisionType PRECLUSTER_THRESHOLD = util::GetUndefined< t_PrecisionType>(),
        const t_PrecisionType SIMILARITY_THRESHOLD = util::GetUndefined< t_PrecisionType>()
      ) :
        Dendrogram< t_DataType, t_PrecisionType>( DISTANCE, DATA, BINARY_PREDICATE, PRECLUSTER_THRESHOLD, SIMILARITY_THRESHOLD)
      {
        util::Stopwatch clustering_timer( " dendrogram clustering ");

        Dendrogram< t_DataType, t_PrecisionType>::m_Node = operator()();
      }

      //! @brief constructor taking ShPtr to t_DataType and a ShPtr to DistanceInterface and a t_PrecisionType
      //! @param DISTANCE is a ShPtr to a distance interface which will define "m_Linkage"
      //! @param DATA is the ShPtrList to t_DataType which will be clustered
      //! @param HEIGHT_CUTOFF is the threshold at which the clustering will be terminated
      //! @param BINARY_PREDICATE defines how two t_DataType are compared
      DendrogramUncached< t_DataType, t_PrecisionType>
      (
        const util::ShPtr< LinkageInterface< t_DataType, t_PrecisionType> > &DISTANCE,
        const util::ShPtr< storage::List< t_DataType> > &DATA,
        const t_PrecisionType HEIGHT_CUTOFF,
        const util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> &BINARY_PREDICATE =
          ( **math::Comparisons< t_PrecisionType>::GetEnums().e_Less),
        const t_PrecisionType PRECLUSTER_THRESHOLD = util::GetUndefined< t_PrecisionType>(),
        const t_PrecisionType SIMILARITY_THRESHOLD = util::GetUndefined< t_PrecisionType>()
      ) :
        Dendrogram< t_DataType, t_PrecisionType>( DISTANCE, DATA, HEIGHT_CUTOFF, BINARY_PREDICATE, PRECLUSTER_THRESHOLD, SIMILARITY_THRESHOLD)
      {
        util::Stopwatch clustering_timer
        (
          " dendrogram clustering ", util::Time( std::numeric_limits< size_t>::max(), 0), util::Message::e_Critical
        );
        Dendrogram< t_DataType, t_PrecisionType>::m_Node = operator()( HEIGHT_CUTOFF);
      }

      //! @brief Clone function
      //! @return pointer to new DendrogramUncached
      DendrogramUncached *Clone() const
      {
        return new DendrogramUncached( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() creates the hierarchy of Nodes
      //! @return returns Node which has the results of the clustering
      Node< t_DataType, t_PrecisionType> operator()()
      {
        //! create Node "node" with the initial base Nodes of clustering
        Node< t_DataType, t_PrecisionType> node( Dendrogram< t_DataType, t_PrecisionType>::InitializeNode()); //< initialize with Nodes made by "m_Data"

        // go until clustering has been completed
        while( node.GetNodes().GetSize() > 2)
        {
          // set "level" to the height at which the next step in clustering occurs
          NextLevel( node);
        }

        // calculate and set the girth of "node" which is the distance between the two nodes it contains
        node.SetGirth
        (
          Dendrogram< t_DataType, t_PrecisionType>::m_Linkage->operator()( **node.GetNodes().Begin(), **( ++node.GetNodes().Begin()))
        );

        // return "node_list" which is the list of Nodes whose girth is less than "HEIGHT_CUTOFF"
        return node;
      }

      //! @brief operator() creates the hierarchy of Nodes
      //! @param HEIGHT_CUTOFF determines the difference between two clusters that stops the agglomeration process
      //!        (inclusive for the cutoff)
      //! @return returns Node which has the results of the clustering
      Node< t_DataType, t_PrecisionType> operator()( const double HEIGHT_CUTOFF)
      {
        //! create Node "node" with the initial base Nodes of clustering
        Node< t_DataType, t_PrecisionType> node( Dendrogram< t_DataType, t_PrecisionType>::InitializeNode()); //< initialize with Nodes made by "m_Data"

        // go until the desired distance/height cutoff is reached
        while
        (
          ( Dendrogram< t_DataType, t_PrecisionType>::m_BinaryPredicate->operator()( node.GetGirth(), HEIGHT_CUTOFF) || !util::IsDefined( node.GetGirth()))
          && node.GetNodes().GetSize() > 2
        )
        {
          // set "level" to the height at which the next step in clustering occurs
          NextLevel( node);
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
        Dendrogram< t_DataType, t_PrecisionType>::Read( ISTREAM);

        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        Dendrogram< t_DataType, t_PrecisionType>::Write( OSTREAM, INDENT);

        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief NextLevel generates the next step in the clustering hierarchy
      //! @param NODE is the list of Nodes which is being clustered together
      //! @return returns a double which is the height of the current level/distance between newly clustered Nodes
      void NextLevel
      (
        Node< t_DataType, t_PrecisionType> &NODE
      )
      {
        // create ShPtrList const_iterators "node_a_itr" and "node_b_itr"
        // these will point to the two Nodes which are nearest to one another so that they can be agglomerated
        typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator node_a_itr, node_b_itr;

        // create double "minimum_linkage" which will hold the minimium distance between any two Nodes
        double minimum_linkage( util::GetUndefined< double>());

        // iterate through "NODES" in order to find the two Nodes which are nearest to one another
        // these two nodes will be agglomerated together to form one node
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
          // increment "itr_b"
          ++itr_b;

          // second iteration through "NODES"
          for
          (
            ;
            itr_b != itr_end;
            ++itr_b
          )
          {
            // create double "current_linkage" which will hold the linkage between the nodes denoted by "itr_a" and "itr_b"
            double current_linkage( Dendrogram< t_DataType, t_PrecisionType>::m_Linkage->operator()( **itr_a, **itr_b));

            // true of Nodes denoted by "itr_a" and "itr_b" are closer than the nodes denoted by
            // "node_a_itr" and "node_b_itr"
            if( !util::IsDefined( minimum_linkage))
            {
              minimum_linkage = current_linkage;
              node_a_itr = itr_a;
              node_b_itr = itr_b;
            }
            else if( Dendrogram< t_DataType, t_PrecisionType>::m_BinaryPredicate->operator()( current_linkage, minimum_linkage))
            {
              minimum_linkage = current_linkage;
              node_a_itr = itr_a;
              node_b_itr = itr_b;
            }
          }
        }

        // "node_a_itr" and "node_b_itr" point to the two nearest Nodes
        // create ShPtrList "combine_nodes"
        util::ShPtrList< Node< t_DataType, t_PrecisionType> > combine_nodes;

        // add Nodes pointed to by "itr_a" and "itr_b" into "combine_nodes"
        combine_nodes.PushBack( *node_a_itr);
        combine_nodes.PushBack( *node_b_itr);

        // create a new Node out of "combine_nodes" with girth of "minimum_linkage" and add this new node to "nodes"
        NODE.GetNodes().PushBack
        (
          util::ShPtr< Node< t_DataType, t_PrecisionType> >( new Node< t_DataType, t_PrecisionType>( combine_nodes, minimum_linkage))
        );

        // remove the Nodes from "nodes" denoted by "itr_a" and "itr_b"
        NODE.GetNodes().Remove( node_a_itr);
        NODE.GetNodes().Remove( node_b_itr);

        // set the girth of "NODE"
        NODE.SetGirth( minimum_linkage);
      }

    private:

    }; // class DendrogramUncached

    // instantiate s_Instance
    template< typename t_DataType, typename t_PrecisionType>
    const util::SiPtr< const util::ObjectInterface> DendrogramUncached< t_DataType, t_PrecisionType>::s_Instance
    (
      GetObjectInstances().AddInstance( new DendrogramUncached< t_DataType, t_PrecisionType>())
    );

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_DENDROGRAM_UNCACHED_H_ 
