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

#ifndef BCL_CLUSTER_OUTPUT_CENTERS_H_
#define BCL_CLUSTER_OUTPUT_CENTERS_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_cluster_node.h"
#include "bcl_cluster_output_rows.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_binary_function_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OutputCenters
    //! @brief TODO: add an general comment to this class
    //!
    //! @see @link example_cluster_output_centers.cpp @endlink
    //! @author alexanns
    //! @date Sep 24, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class OutputCenters :
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
      OutputCenters() :
        m_MemberDistanceFunction(),
        m_BinaryPredicate()
      {
      }

      //! @brief constructor taking all member variables as parameters
      //! @param MEMBER_DISTANCE_FUNCTION
      //! @param BINARY_PREDICATE
      OutputCenters
      (
        const util::ShPtr< math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const t_DataType> >, t_PrecisionType> >
          &MEMBER_DISTANCE_FUNCTION,
        const util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > &BINARY_PREDICATE
      ) :
        m_MemberDistanceFunction( MEMBER_DISTANCE_FUNCTION),
        m_BinaryPredicate( BINARY_PREDICATE)
      {
      }

      //! @brief Clone function
      //! @return pointer to new OutputCenters
      OutputCenters< t_DataType, t_PrecisionType> *Clone() const
      {
        return new OutputCenters< t_DataType, t_PrecisionType>( *this);
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
        io::Serialize::Write( m_BinaryPredicate, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief OutputNodeMembers outputs the individual members of a Node
      //! @param OFSTREAM the output stream to which the Node will be written
      //! @param NODE the node which will be output to OFSTREAM
      //! @param NODE_COUNTER the number of the current node
      //! @param OUTER_NODE_COUNTER
      //! @return returns io::OFStream
      std::ostream &OutputNodeMembers
      (
        std::ostream &OFSTREAM, const Node< t_DataType, t_PrecisionType> &NODE
      ) const
      {
        // create t_PrecisionType "girth" and initialize with the girth of "NODE"
        const t_PrecisionType girth( NODE.GetGirth());

        // create const size_t "size" and initialize with the number of members in "NODE"
        const size_t size( NODE.GetMembers().GetSize());

        // create SiPtr "node_center" to the center of "NODE"
        const util::SiPtr< const t_DataType> node_center
        (
          NODE.GetCenter( m_MemberDistanceFunction, m_BinaryPredicate)
        );

        bool leaf( NODE.IsLeaf());

        // output the node number and the girth of the current node denoted by "itr"
        OFSTREAM << "NODE ";
        // write out the identifier; it might be nan, so use io::Serialize
        io::Serialize::Write( NODE.GetIdentifier(), OFSTREAM, 0);
        OFSTREAM << " : Member : " << std::string( util::Format()( *node_center))
                 << " : Size : " << size
                 << " : Leaf : " << util::Format()( leaf)
                 << " : Linkage : ";
        // write out the girth; it might be nan, so use io::Serialize
        io::Serialize::Write( girth, OFSTREAM, 0) << '\n';

        // end
        return OFSTREAM;
      }

    }; // class OutputCenters

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_OUTPUT_CENTERS_H_ 
