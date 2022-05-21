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

#ifndef BCL_CHEMISTRY_CONFORMATION_SET_SAME_CONFIGURATION_H_
#define BCL_CHEMISTRY_CONFORMATION_SET_SAME_CONFIGURATION_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_comparison_interface.h"
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_fragment_configuration_shared.h"
#include "bcl_chemistry_fragment_conformation_shared.h"
#include "graph/bcl_graph_const_graph.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConformationSetSameConfiguration
    //! @brief Container class for FragmentConstitutionShared objects
    //!
    //! @see @link example_chemistry_conformation_set_same_configuration.cpp @endlink
    //! @author kothiwsk
    //! @date Mar 12, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConformationSetSameConfiguration :
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
      double                                                 m_Tolerance;

      //! graph of constitution
      graph::ConstGraph< size_t, size_t>                     m_ConfigurationGraph;

      //! list of FragmentConstitutionShared
      util::ShPtr< FragmentConfigurationShared>              m_Configuration;

      //! list of FragmentConfigurationShared
      util::ShPtrList< FragmentConformationShared>           m_Conformations;

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
      ConformationSetSameConfiguration();

      //! @brief constructor
      //! @param METHOD conformation comparison method to determine if conformations are same
      //! @param TOLERANCE tolerance value for comparing conformation values
      ConformationSetSameConfiguration
      (
        util::Implementation< ConformationComparisonInterface> METHOD,
        double TOLERANCE
      );

      //! @brief Clone function
      //! @return pointer to new ConformationSetSameConfiguration
      ConformationSetSameConfiguration *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns the configuration of the set
      //! @return configuration of the set
      const util::ShPtr< FragmentConfigurationShared> &GetConfiguration() const
      {
        return m_Configuration;
      }

      //! @brief return a list containing the conformations
      //! @return list of conformations
      const util::ShPtrList< FragmentConformationShared> &GetConformations() const
      {
        return m_Conformations;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns the number of conformations
      //! @return the number of conformations
      size_t GetSize() const
      {
        return m_Conformations.GetSize();
      }

      //! @brief const iterator begin
      //! @return iterator begin of container storing conformations
      const_iterator Begin() const
      {
        return m_Conformations.Begin();
      }

      //! @brief const iterator end
      //! @return iterator end of container storing conformations
      const_iterator End() const
      {
        return m_Conformations.End();
      }

      //! @brief append a molecule to configuration set
      //! @param FRAGMENT fragment conformation shared that needs to be added to configuration set
      //! @param CONFIGURATION configuration of interest
      //! @param ISOMORPHISM isomorphism which if it was determined at constitution layer
      //! @return a pair of iterator to conformation and a bool, true indicating conformation insert was successful
      std::pair< const_iterator, bool> Insert
      (
        const ConformationInterface &FRAGMENT,
        util::ShPtr< FragmentConfigurationShared> CONFIGURATION = util::ShPtr< FragmentConfigurationShared>(),
        util::SiPtr< storage::Vector< size_t> > ISOMORPHISM = util::SiPtr< storage::Vector< size_t> >()
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

    }; // class ConformationSetSameConfiguration

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFORMATION_SET_SAME_CONFIGURATION_H_
