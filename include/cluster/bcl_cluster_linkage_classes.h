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

#ifndef BCL_CLUSTER_LINKAGE_CLASSES_H_
#define BCL_CLUSTER_LINKAGE_CLASSES_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_cluster_linkage_average.h"
#include "bcl_cluster_linkage_complete.h"
#include "bcl_cluster_linkage_single.h"
#include "bcl_cluster_linkage_total.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LinkageClasses
    //! @brief TODO: document
    //!
    //! @see @link example_cluster_linkage_classes.cpp @endlink
    //! @author alexanns
    //! @date 10/01/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class LinkageClasses :
      public util::Enumerate< util::ShPtr< LinkageInterface< t_DataType, t_PrecisionType> >, LinkageClasses< t_DataType, t_PrecisionType> >
    {
      friend class util::Enumerate< util::ShPtr< LinkageInterface< t_DataType, t_PrecisionType> >, LinkageClasses< t_DataType, t_PrecisionType> >;

    public:

      // typedef for LinkageClass
      typedef typename util::Enumerate< util::ShPtr< LinkageInterface< t_DataType, t_PrecisionType> >, LinkageClasses< t_DataType, t_PrecisionType> >::EnumType LinkageClass;
      typedef typename util::Enumerate< util::ShPtr< LinkageInterface< t_DataType, t_PrecisionType> >, LinkageClasses< t_DataType, t_PrecisionType> >::EnumDataType LinkageClassData;

    //////////
    // data //
    //////////

      // declare all linkage classes //
      const LinkageClass e_Complete;
      const LinkageClass e_Single;
      const LinkageClass e_Average;
      const LinkageClass e_Total;

      using util::Enumerate< util::ShPtr< LinkageInterface< t_DataType, t_PrecisionType> >, LinkageClasses< t_DataType, t_PrecisionType> >::AddEnum;

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct all LinkageClasses
      LinkageClasses< t_DataType, t_PrecisionType>() :
        e_Complete( AddEnum( "Complete", util::ShPtr< LinkageInterface< t_DataType, t_PrecisionType> >( new LinkageComplete< t_DataType, t_PrecisionType>()))),
        e_Single(   AddEnum( "Single",   util::ShPtr< LinkageInterface< t_DataType, t_PrecisionType> >( new LinkageSingle< t_DataType, t_PrecisionType>()))),
        e_Average(  AddEnum( "Average",  util::ShPtr< LinkageInterface< t_DataType, t_PrecisionType> >( new LinkageAverage< t_DataType, t_PrecisionType>()))),
        e_Total(    AddEnum( "Total",    util::ShPtr< LinkageInterface< t_DataType, t_PrecisionType> >( new LinkageTotal< t_DataType, t_PrecisionType>())))
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

    }; // class LinkageClasses

    //! @brief construct on access function for all LinkageClasses
    //! @return reference to only instances of LinkageClasses
    template< typename t_DataType, typename t_PrecisionType> inline const LinkageClasses< t_DataType, t_PrecisionType> &GetLinkageClasses()
    {
      return LinkageClasses< t_DataType, t_PrecisionType>::GetEnums();
    }

  } // namespace cluster

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< cluster::LinkageInterface< std::string, float> >, cluster::LinkageClasses< std::string, float> >;

  } // namespace util
} // namespace bcl

#endif //BCL_CLUSTER_LINKAGE_CLASSES_H_
