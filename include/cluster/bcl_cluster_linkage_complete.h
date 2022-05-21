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

#ifndef BCL_CLUSTER_LINKAGE_COMPLETE_H_
#define BCL_CLUSTER_LINKAGE_COMPLETE_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_cluster_linkage_interface.h"
#include "bcl_cluster_node.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LinkageComplete
    //! @brief LinkageComplete determines the distance between two clusters as the two points that are farthest between
    //!        the two clusters.
    //!
    //! @see @link example_cluster_linkage_complete.cpp @endlink
    //! @author alexanns, woetzen
    //! @date June 6, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class LinkageComplete :
      public LinkageInterface< t_DataType, t_PrecisionType>
    {
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
      LinkageComplete() :
        LinkageInterface< t_DataType, t_PrecisionType>()
      {
      }

      //! @brief constructor from a FunctionInterface
      LinkageComplete
      (
        const util::ShPtr< math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const t_DataType> >, t_PrecisionType> > &FUNCTION,
        const util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > &BINARY_PREDICATE
      ) :
        LinkageInterface< t_DataType, t_PrecisionType>( FUNCTION, BINARY_PREDICATE)
      {
      }

      //! @brief copy constructor
      //! @param LINKAGE_COMPLETE the LinkageComplete< t_DataType> from which this LinkageComplete will be copied
      LinkageComplete( const LinkageComplete< t_DataType, t_PrecisionType> &LINKAGE_COMPLETE) :
        LinkageInterface< t_DataType, t_PrecisionType>( LINKAGE_COMPLETE)
      {
      }

      //! @brief Clone function
      //! @return pointer to new LinkageComplete< t_DataType>
      LinkageComplete< t_DataType, t_PrecisionType> *Clone() const
      {
        return new LinkageComplete< t_DataType, t_PrecisionType>( *this);
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

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() which takes two Nodes and returns the distance between them as the two farthest points
      //! @param NODE_A is the first Node involved in the distance calculation
      //! @param NODE_B is the second Node involved in the distance calculation
      //! @return returns a t_PrecisionType which is the distance between the two clusters
      t_PrecisionType operator()( const Node< t_DataType, t_PrecisionType> &NODE_A, const Node< t_DataType, t_PrecisionType> &NODE_B) const
      {
        // create t_PrecisionType "max_distance" to hold the furthest distance between any two members of "NODE_A" and "NODE_B"
        t_PrecisionType max_distance( util::GetUndefined< t_PrecisionType>()); //< initialize to zero

        // iterate through "NODE_A"
        for
        (
          typename util::SiPtrList< t_DataType>::const_iterator
           node_a_itr( NODE_A.GetMembers().Begin()), node_a_itr_end( NODE_A.GetMembers().End());
          node_a_itr != node_a_itr_end;
          ++node_a_itr
        )
        {
          // iterate through "NODE_B"
          for
          (
            typename util::SiPtrList< t_DataType>::const_iterator
              node_b_itr( NODE_B.GetMembers().Begin()), node_b_itr_end( NODE_B.GetMembers().End());
            node_b_itr != node_b_itr_end;
            ++node_b_itr
          )
          {
            // create t_PrecisionType "current_distance" and initialize to the distance between the t_DataTypes pointed
            // by to "node_a_itr" and "node_b_itr"
            const t_PrecisionType current_distance
            (
              LinkageInterface< t_DataType, t_PrecisionType>::m_Function->operator()( storage::VectorND< 2, util::SiPtr< const t_DataType> >( **node_a_itr, **node_b_itr))
            );

            // true if "max_distance" is not defined - means can just assign it to whatever "current_distance" is
            if( !util::IsDefined( max_distance))
            {
              max_distance = current_distance;
            }
            // true if "current_distance" is larger than "max_distance"
            else if( !LinkageInterface< t_DataType, t_PrecisionType>::m_BinaryPredicate->operator()( current_distance, max_distance))
            {
              max_distance = current_distance;
            }
          }
        }

        // return "max_distance" which is the distance between "NODE_A" and "NODE_B"
        return max_distance;
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
        // read "m_Function"
        io::Serialize::Read( LinkageInterface< t_DataType, t_PrecisionType>::m_Function, ISTREAM);

        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write "m_Function"
        io::Serialize::Write( LinkageInterface< t_DataType, t_PrecisionType>::m_Function, OSTREAM, INDENT);

        return OSTREAM;
      }

    }; // LinkageComplete

    // instantiate s_Instance
    template< typename t_DataType, typename t_PrecisionType>
    const util::SiPtr< const util::ObjectInterface> LinkageComplete< t_DataType, t_PrecisionType>::s_Instance
    (
      GetObjectInstances().AddInstance( new LinkageComplete< t_DataType, t_PrecisionType>())
    );

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_LINKAGE_COMPLETE_H_ 
