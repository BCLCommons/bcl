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

#ifndef BCL_CLUSTER_OUTPUT_SORTED_MATRIX_H_
#define BCL_CLUSTER_OUTPUT_SORTED_MATRIX_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_cluster_node.h"
#include "bcl_cluster_output_rows.h"
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_binary_function_interface.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OutputSortedMatrix
    //! @brief TODO: add an general comment to this class
    //!
    //! @see @link example_cluster_output_sorted_matrix.cpp @endlink
    //! @author mendenjl
    //! @date Oct 29, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class OutputSortedMatrix :
      public OutputRows< t_DataType, t_PrecisionType>
    {

    private:

    //////////
    // data //
    //////////

      //! "m_MemberDistanceFunction" is the function which is used to compare the individual members of the Nodes
      util::ShPtr
      <
        math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const t_DataType> >, t_PrecisionType>
      > m_MemberDistanceFunction;

      //! ShPtr to a BinaryFunctionInterface "m_BinaryPredicate" to defines how two t_DataType are compared
      util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > m_BinaryPredicate;

      //! bool, controls whether to output as a table or a matrix
      bool m_TabularOutput;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Default constructor
      OutputSortedMatrix( const bool &TABULAR_OUTPUT = false) :
        m_MemberDistanceFunction(),
        m_BinaryPredicate(),
        m_TabularOutput( TABULAR_OUTPUT)
      {
      }

      //! @brief constructor taking all member variables as parameters
      //! @param MEMBER_DISTANCE_FUNCTION
      //! @param BINARY_PREDICATE
      OutputSortedMatrix
      (
        const util::ShPtr< math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const t_DataType> >, t_PrecisionType> >
          &MEMBER_DISTANCE_FUNCTION,
        const util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > &BINARY_PREDICATE,
        const bool &TABULAR_OUTPUT
      ) :
        m_MemberDistanceFunction( MEMBER_DISTANCE_FUNCTION),
        m_BinaryPredicate( BINARY_PREDICATE),
        m_TabularOutput( TABULAR_OUTPUT)
      {
      }

      //! @brief Clone function
      //! @return pointer to new OutputSortedMatrix
      OutputSortedMatrix< t_DataType, t_PrecisionType> *Clone() const
      {
        return new OutputSortedMatrix< t_DataType, t_PrecisionType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief WriteOutput writes some data of a list of Nodes to some output stream
      //! @param FILENAME is the stream to which the output will be written
      //! @param NODE_LIST the list of nodes which is going to be output
      //! @return returns data construct holding the distances between all objects for use in a LinkageInterface
      void WriteOutput
      (
        const std::string &FILENAME, const util::SiPtrList< const Node< t_DataType, t_PrecisionType> > &NODE_LIST
      ) const
      {
        io::OFStream ofstream;
        io::File::MustOpenOFStream( ofstream, FILENAME);
        util::SiPtrVector< const t_DataType> node_members_sorted;
        {
          // get the expanded node
          util::SiPtrList< const Node< t_DataType, t_PrecisionType> > internal_nodes
          (
            NODE_LIST.FirstElement()->ExpandAllNodes()
          );

          // create a list with all the members, which will now be ordered by their clustering position
          util::SiPtrList< const t_DataType> node_members_sorted_list;

          // get the members, ordered by the node ordering
          for
          (
            typename util::SiPtrList< const Node< t_DataType, t_PrecisionType> >::const_iterator
              itr( internal_nodes.Begin()), itr_end( internal_nodes.End());
            itr != itr_end;
            ++itr
          )
          {
            node_members_sorted_list.Append( ( *itr)->GetOrphanedMembers());
          }
          node_members_sorted =
            util::SiPtrVector< const t_DataType>( node_members_sorted_list.Begin(), node_members_sorted_list.End());
        }

        {
          // write the table headings, if desired
          if( m_TabularOutput)
          {
            ofstream << "Table\t";
            for
            (
              typename util::SiPtrVector< const t_DataType>::const_iterator
                itr_a( node_members_sorted.Begin()), itr_end( node_members_sorted.End());
              itr_a != itr_end;
              ++itr_a
            )
            {
              ofstream << **itr_a << '\t';
            }
            ofstream << '\n';
          }

          for
          (
            typename util::SiPtrVector< const t_DataType>::const_iterator
              itr_a( node_members_sorted.Begin()), itr_end( node_members_sorted.End());
            itr_a != itr_end;
            ++itr_a
          )
          {
            if( m_TabularOutput)
            {
              ofstream << **itr_a << '\t';
            }
            for
            (
              typename util::SiPtrVector< const t_DataType>::const_iterator itr_b( node_members_sorted.Begin());
              itr_b != itr_end;
              ++itr_b
            )
            {
              if( itr_b != itr_a)
              {
                const t_PrecisionType distance_similarity
                (
                  m_MemberDistanceFunction->operator()
                  (
                    storage::VectorND< 2, util::SiPtr< const t_DataType> >( *itr_a, *itr_b)
                  )
                );
                ofstream << distance_similarity << '\t';
              }
              else if( m_BinaryPredicate == *math::Comparisons< t_PrecisionType>::GetEnums().e_Less)
              {
                ofstream << "0\t";
              }
              else
              {
                ofstream << "1\t";
              }
            }
            ofstream << '\n';
          }
        }
        io::File::CloseClearFStream( ofstream);
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
        io::Serialize::Read( m_MemberDistanceFunction, ISTREAM);
        io::Serialize::Read( m_BinaryPredicate, ISTREAM);
        io::Serialize::Read( m_TabularOutput, ISTREAM);

        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_MemberDistanceFunction, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_BinaryPredicate, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_TabularOutput, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

    }; // class OutputSortedMatrix

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_OUTPUT_SORTED_MATRIX_H_
