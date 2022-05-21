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

#ifndef BCL_CLUSTER_LINKAGE_INTERFACE_H_
#define BCL_CLUSTER_LINKAGE_INTERFACE_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
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
    //! @class LinkageInterface
    //! @brief LinkageInterface is the interface class from which methods for calculating the linkage between Nodes
    //!        should be derived.
    //!
    //! @remarks example unnecessary
    //! @author alexanns, woetzen
    //! @date June 6, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class LinkageInterface :
      public util::ObjectInterface
    {
    protected:

    //////////
    // data //
    //////////

      //! "m_function" is the function which is used to compare the individual members of the Nodes
      util::ShPtr
      <
        math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const t_DataType> >, t_PrecisionType>
      > m_Function;

      //! ShPtr to a BinaryFunctionInterface "m_BinaryPredicate" to defines how linkages are determined
      util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > m_BinaryPredicate;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LinkageInterface() :
        m_Function(),
        m_BinaryPredicate()
      {
      }

      //! @brief constructor from a FunctionInterface and a BinaryFunctionInterface
      //! @param FUNCTION the function which is used to compare the individual members of the Nodes
      LinkageInterface
      (
        const util::ShPtr
        <
          math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const t_DataType> >, t_PrecisionType>
        > &FUNCTION
      ) :
        m_Function( FUNCTION),
        m_BinaryPredicate()
      {
      }

      //! @brief constructor from a FunctionInterface and a BinaryFunctionInterface
      //! @param FUNCTION the function which is used to compare the individual members of the Nodes
      //! @param BINARY_PREDICATE
      LinkageInterface
      (
        const util::ShPtr
        <
          math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const t_DataType> >, t_PrecisionType>
        > &FUNCTION,
        const util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > &BINARY_PREDICATE
      ) :
        m_Function( FUNCTION),
        m_BinaryPredicate( BINARY_PREDICATE)
      {
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief SetDistanceFunction sets "DISTANCE_FUNCTION" to some distance function
      //! @param DISTANCE_FUNCTION the distance fucntion which "m_Function" will be set to
      void SetDistanceFunction
      (
        const util::ShPtr< math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const t_DataType> >, t_PrecisionType> >
        &DISTANCE_FUNCTION
      )
      {
        m_Function = DISTANCE_FUNCTION;
      }

      //! @brief GetDistanceFunction returns a const reference to "m_Function"
      //! @return returns a const reference to "m_Function"
      const util::ShPtr< math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const t_DataType> >, t_PrecisionType> > &
      GetDistanceFunction() const
      {
        return m_Function;
      }

      void SetBinaryPredicate( const util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> &BINARY_PREDICATE)
      {
        m_BinaryPredicate =
          util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> >( BINARY_PREDICATE.Clone());
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() which takes two Nodes and returns the linkage between them
      //! @param NODE_A is the first Node involved in the linkage calculation
      //! @param NODE_B is the second Node involved in the linkage calculation
      //! @return returns a t_PrecisionType which is the linkage between the two clusters
      virtual t_PrecisionType operator()( const Node< t_DataType, t_PrecisionType> &NODE_A, const Node< t_DataType, t_PrecisionType> &NODE_B) const = 0;

      //! @brief operator() which takes a Node and returns the linkage within that existing cluster
      //! @param NODE is the Node to calculate the linkage for
      //! @return returns a t_PrecisionType which is the linkage in this clusters
      virtual t_PrecisionType operator()( const Node< t_DataType, t_PrecisionType> &NODE) const
      {
        return this->operator()( NODE, NODE);
      }

    }; // LinkageInterface

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_LINKAGE_INTERFACE_H_
