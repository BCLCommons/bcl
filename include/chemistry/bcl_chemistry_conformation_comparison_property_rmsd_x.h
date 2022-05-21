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

#ifndef BCL_CHEMISTRY_CONFORMATION_COMPARISON_PROPERTY_RMSD_X_H_
#define BCL_CHEMISTRY_CONFORMATION_COMPARISON_PROPERTY_RMSD_X_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_comparison_interface.h"
#include "descriptor/bcl_descriptor_base.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConformationComparisonPropertyRMSDX
    //! @brief This class is designed to be used for determining and comparing 3D structures for molecules based
    //!        on euclidean and property distance, and even for molecules with differing constitutions
    //!
    //! @see @link example_chemistry_conformation_comparison_property_rmsd_x.cpp @endlink
    //! @author mendenjl
    //! @date Apr 12, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConformationComparisonPropertyRMSDX :
      public ConformationComparisonInterface
    {

    protected:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! Atom properties to include in the distance
      mutable storage::Vector< descriptor::CheminfoProperty> m_Properties;

      //! weight to give the atom property
      double m_PropertyWeight;

      //! Base - approximate # of atoms that should go into an alignment
      double m_Base;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ConformationComparisonPropertyRMSDX();

      //! virtual copy constructor
      ConformationComparisonPropertyRMSDX *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief find the RMSD between two conformations
      //! @param MOLECULE_B - the fragment to align against
      //! @param MOLECULE_B - the fragment being aligned
      //! @return the RMSD between MOLECULE_A and MOLECULE_B
      double operator()
      (
        const ConformationInterface &MOLECULE_A,
        const ConformationInterface &MOLECULE_B
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief prepare the class for comparing a conformation
      //! @param MOLECULE the molecule to prepare to compare
      void Prepare( const ConformationInterface &MOLECULE) const;

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    };

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFORMATION_COMPARISON_PROPERTY_RMSD_X_H_
