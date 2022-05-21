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

#ifndef BCL_APP_PHARM_MAP_H
#define BCL_APP_PHARM_MAP_H
// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "app/bcl_app.h"
#include "app/bcl_app_interface.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_interface.h"
#include "graph/bcl_graph_common_subgraph_isomorphism.h"
#include "linal/bcl_linal_vector.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"

// include headers from the bcl - sorted alphabetically

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PharmMap
    //! @brief Application for generating libraries for synthesis using QSAR models and random structure generator
    //!
    //! @author geanesar
    //! @date 08/31/2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class PharmMap :
      public Interface
    {
    private:

    //////////
    // data //
    //////////

      //! flag sets which method, building or scoring, the method will do
      util::ShPtr< command::FlagInterface> m_MethodFlag;

      //! flag sets the filename for the scaffold
      util::ShPtr< command::FlagInterface> m_ScaffoldFlag;

      //! flag sets where derivatives of base fragment are stored
      util::ShPtr< command::FlagInterface> m_FragmentsFlag;

      //! the number of derivatives to generate
      util::ShPtr< command::FlagInterface> m_SampleSizeFlag;

      //! flag sets where to get a set of generated molecules from
      util::ShPtr< command::FlagInterface> m_DerivativesFlag;

      //! flag sets the output sdf filename
      util::ShPtr< command::FlagInterface> m_PymolOutputFilenameFlag;

      //! flag sets the output sdf filename
      util::ShPtr< command::FlagInterface> m_CSVOutputFilenameFlag;

      //! flag sets the filename containing object data labels of properties to feed to PharmMap
      util::ShPtr< command::FlagInterface> m_PropertiesToMapFlag;

      //! flag that takes the models to use
      util::ShPtr< command::FlagInterface> m_ScorerFlag;

      //! flag sets where to output generated molecules, if any
      util::ShPtr< command::FlagInterface> m_MoleculeOutputFlag;

      //! how many grow points to change
      util::ShPtr< command::FlagInterface> m_RemoveZerodPropertiesFlag;

      //! where to write detailed output to
      util::ShPtr< command::FlagInterface> m_DetailsFilenameFlag;

      //! file where pharmmap coefficients are stored
      util::ShPtr< command::FlagInterface> m_PharmMapFileFlag;

      //! whether score pairs are treated in a binary fashion
      util::ShPtr< command::FlagInterface> m_BinaryDifferenceFlag;

      //! score tolerance
      util::ShPtr< command::FlagInterface> m_ScoreToleranceFlag;

      util::ShPtr< command::FlagInterface> m_PercentCoverageForScaffoldFlag;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! @brief Default constructor
      PharmMap();

    public:

      //! @brief Clone function
      //! @return pointer to new PharmMap
      PharmMap *Clone() const
      {
        return new PharmMap( *this);
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

      //! @brief initializes the command object for this application
      //! @return a ShPtr to a Command containing all of this applications parameters
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return 0 for success
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief a helper class for parsing a CSV file containing pharmmap info
      storage::Map< size_t, linal::Vector< float> > ParsePharmMapCSV( const std::string &FILENAME) const;

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      // Static instance of PharmMap
      static const ApplicationType PharmMap_Instance;

    }; // class PharmMap

  } // namespace app
} // namespace bcl

#endif // BCL_APP_PHARM_MAP_H
