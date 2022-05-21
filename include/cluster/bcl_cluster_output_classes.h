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

#ifndef BCL_CLUSTER_OUTPUT_CLASSES_H_
#define BCL_CLUSTER_OUTPUT_CLASSES_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_cluster_output_centers.h"
#include "bcl_cluster_output_sorted_matrix.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OutputClasses
    //! @brief OutputClasses is an enumerator for the formats/methods for outputting the results of clustering.
    //!
    //! @see @link example_cluster_output_classes.cpp @endlink
    //! @author alexanns
    //! @date October 09, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class OutputClasses :
      public util::Enumerate< util::ShPtr< OutputInterface< t_DataType, t_PrecisionType> >, OutputClasses< t_DataType, t_PrecisionType> >
    {
      friend class util::Enumerate< util::ShPtr< OutputInterface< t_DataType, t_PrecisionType> >, OutputClasses< t_DataType, t_PrecisionType> >;

    public:

      // typedef for OutputClass
      typedef typename util::Enumerate< util::ShPtr< OutputInterface< t_DataType, t_PrecisionType> >, OutputClasses>::EnumType OutputClass;
      typedef typename util::Enumerate< util::ShPtr< OutputInterface< t_DataType, t_PrecisionType> >, OutputClasses>::EnumDataType OutputClassData;

    //////////
    // data //
    //////////

      // declare all input classes //

      //! bcl::storage::Table formatted input; all information in table is needed
      const OutputClass e_Rows;
      const OutputClass e_Centers;
      const OutputClass e_Matrix;
      const OutputClass e_Table;

      using util::Enumerate< util::ShPtr< OutputInterface< t_DataType, t_PrecisionType> >, OutputClasses< t_DataType, t_PrecisionType> >::AddEnum;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all InputClasses
      OutputClasses() :
        e_Rows( AddEnum( "Rows", util::ShPtr< OutputInterface< t_DataType, t_PrecisionType> >( new OutputRows< t_DataType, t_PrecisionType>()))),
        e_Centers( AddEnum( "Centers", util::ShPtr< OutputInterface< t_DataType, t_PrecisionType> >( new OutputCenters< t_DataType, t_PrecisionType>()))),
        e_Matrix( AddEnum( "Matrix", util::ShPtr< OutputInterface< t_DataType, t_PrecisionType> >( new OutputSortedMatrix< t_DataType, t_PrecisionType>( false)))),
        e_Table( AddEnum( "Table", util::ShPtr< OutputInterface< t_DataType, t_PrecisionType> >( new OutputSortedMatrix< t_DataType, t_PrecisionType>( true))))
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

    }; // class OutputClasses

    //! @brief construct on access function for all InputClasses
    //! @return reference to only instances of InputClasses
    template< typename t_DataType, typename t_PrecisionType>
    inline const OutputClasses< t_DataType, t_PrecisionType> &GetOutputClasses()
    {
      return OutputClasses< t_DataType, t_PrecisionType>::GetEnums();
    }

  } // namespace cluster

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< cluster::OutputInterface< std::string, float> >, cluster::OutputClasses< std::string, float> >;
    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< cluster::OutputInterface< linal::Vector< float>, float> >, cluster::OutputClasses< linal::Vector< float>, float> >;

  } // namespace util
} // namespace bcl

#endif //BCL_CLUSTER_OUTPUT_CLASSES_H_
