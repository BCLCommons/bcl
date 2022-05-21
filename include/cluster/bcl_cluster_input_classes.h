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

#ifndef BCL_CLUSTER_INPUT_CLASSES_H_
#define BCL_CLUSTER_INPUT_CLASSES_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_cluster_input_pairwise_list.h"
#include "bcl_cluster_input_table.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class InputClasses
    //! @brief for different formats that can be read in for clustering
    //!
    //! @see @link example_cluster_input_classes.cpp @endlink
    //! @author alexanns
    //! @date 10/01/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class InputClasses :
      public util::Enumerate< util::ShPtr< InputInterface< t_DataType, t_PrecisionType> >, InputClasses< t_DataType, t_PrecisionType> >
    {
      friend class util::Enumerate< util::ShPtr< InputInterface< t_DataType, t_PrecisionType> >, InputClasses< t_DataType, t_PrecisionType> >;

    public:

      // typedef for InputClass
      typedef typename util::Enumerate< util::ShPtr< InputInterface< t_DataType, t_PrecisionType> >, InputClasses>::EnumType InputClass;
      typedef typename util::Enumerate< util::ShPtr< InputInterface< t_DataType, t_PrecisionType> >, InputClasses>::EnumDataType InputClassData;

    //////////
    // data //
    //////////

      // declare all input classes

       //! bcl::storage::Table formatted input but only the upper triangle information is needed
      const InputClass e_TableUpperTriangle;
      const InputClass e_TableLowerTriangle;
      const InputClass e_TableComplete;

      const InputClass e_PairwiseList;

      using util::Enumerate< util::ShPtr< InputInterface< t_DataType, t_PrecisionType> >, InputClasses< t_DataType, t_PrecisionType> >::AddEnum;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all InputClasses
      InputClasses() :
        e_TableUpperTriangle( AddEnum( "TableUpperTriangle", util::ShPtr< InputInterface< t_DataType, t_PrecisionType> >( new InputTable< t_PrecisionType>( true, false)))),
        e_TableLowerTriangle( AddEnum( "TableLowerTriangle", util::ShPtr< InputInterface< t_DataType, t_PrecisionType> >( new InputTable< t_PrecisionType>( false, true)))),
        e_TableComplete(      AddEnum( "TableComplete",      util::ShPtr< InputInterface< t_DataType, t_PrecisionType> >( new InputTable< t_PrecisionType>( true, true)))),
        e_PairwiseList(       AddEnum( "PairwiseList",       util::ShPtr< InputInterface< t_DataType, t_PrecisionType> >( new InputPairwiseList< t_PrecisionType>())))
      {
      }

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    }; // class InputClasses

    //! @brief construct on access function for all InputClasses
    //! @return reference to only instances of InputClasses
    template< typename t_DataType, typename t_PrecisionType> inline const InputClasses< t_DataType, t_PrecisionType> &GetInputClasses()
    {
      return InputClasses< t_DataType, t_PrecisionType>::GetEnums();
    }

  } // namespace cluster

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< cluster::InputInterface< std::string, float> >, cluster::InputClasses< std::string, float> >;
    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< cluster::InputInterface< linal::Vector< float>, float> >, cluster::InputClasses< linal::Vector< float>, float> >;

  } // namespace util
} // namespace bcl

#endif //BCL_CLUSTER_INPUT_CLASSES_H_
