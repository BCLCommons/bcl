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

#ifndef BCL_CHEMISTRY_CONFORMATION_SET_SAME_CONSTITUTION_H_
#define BCL_CHEMISTRY_CONFORMATION_SET_SAME_CONSTITUTION_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_configuration_set_same_constitution.h"
#include "bcl_chemistry_conformation_comparison_interface.h"
#include "bcl_chemistry_conformation_set_same_configuration.h"
#include "bcl_chemistry_fragment_conformation_shared.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConformationSetSameConstitution
    //! @brief Container class for FragmentConstitutionShared objects
    //!
    //! @see @link example_chemistry_conformation_set_same_constitution.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Mar 14, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConformationSetSameConstitution :
      public util::ObjectInterface
    {
    public:

    //////////////
    // typedefs //
    //////////////

      typedef util::ShPtrList< FragmentConformationShared>::const_iterator const_iterator;

    private:

    //////////
    // data //
    //////////

      //! metric to compare conformations
      util::Implementation< ConformationComparisonInterface> m_ConfomationComparer;

      //! tolerance value to determine if the conformations are different
      double                                                  m_Tolerance;

      //! All conformations
      util::ShPtrList< FragmentConformationShared>            m_Conformations;

      //! set to store configurations that belong to same constitution
      ConfigurationSetSameConstitution                        m_ConfigurationSet;

      //! map of constitution to various configurations
      storage::Map
      <
        util::SiPtr< const FragmentConfigurationShared>,
        ConformationSetSameConfiguration
      >                                                       m_ConformationsMap;

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
      ConformationSetSameConstitution();

      //! @brief constructor
      //! @param METHOD conformation comparison method to determine if conformations are same
      //! @param TOLERANCE tolerance value for comparing conformation values
      ConformationSetSameConstitution
      (
        util::Implementation< ConformationComparisonInterface> METHOD,
        double TOLERANCE = double( 0.5)
      );

      //! @brief Clone function
      //! @return pointer to new ConformationSetSameConstitution
      ConformationSetSameConstitution *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return a list containing the constitutions
      //! @return list of constitutions
      const ConfigurationSetSameConstitution &GetConfigurations() const
      {
        return m_ConfigurationSet;
      }

      //! @brief return a list containing the conformations
      //! @return list of conformations
      const util::ShPtrList< FragmentConformationShared> &GetConformations() const
      {
        return m_Conformations;
      }

      //! @brief return a reference to ConformationsMap
      //! @return map of ConformationsMap
      const storage::Map< util::SiPtr< const FragmentConfigurationShared>, ConformationSetSameConfiguration> &GetConformationsMap() const
      {
        return m_ConformationsMap;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns the number of conformations
      size_t GetSize() const
      {
        return m_Conformations.GetSize();
      }

      //! @brief const iterator begin
      const_iterator Begin() const
      {
        return m_Conformations.Begin();
      }

      //! @brief const iterator end
      const_iterator End() const
      {
        return m_Conformations.End();
      }

      //! @brief append a molecule to SetConstitution
      //! @param FRAGMENT fragment conformation  that needs to be added to set constitution
      //! @param CONSTITUTION constitution to link the conformations with
      //! @param ISOMORPHISM isomorphism which if it was determined at constitution layer
      //! @return a pair of iterator to conformation and a bool, true indicating conformation insert was successful
      std::pair< const_iterator, bool> Insert
      (
        const ConformationInterface &FRAGMENT,
        const util::ShPtr< FragmentConstitutionShared> &CONSTITUTION,
        storage::Vector< size_t> &ISOMORPHISM
      );

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class ConformationSetSameConstitution
  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFORMATION_SET_SAME_CONSTITUTION_H_
