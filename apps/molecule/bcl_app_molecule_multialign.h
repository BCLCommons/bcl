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

#ifndef BCL_APP_MOLECULE_MULTIALIGN_H_
#define BCL_APP_MOLECULE_MULTIALIGN_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"
#include "graph/bcl_graph.fwd.hh"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_comparison_multi_align.h"
#include "command/bcl_command_command.h"
#include "linal/bcl_linal_matrix.h"
#include "sched/bcl_sched_mutex.h"
#include "storage/bcl_storage_pair.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeMultiAlign
    //! @brief Application for performing simultaneous multiple flexible alignment of > 2 small molecules
    //!
    //! @see @link example_app_molecule_multialign.cpp @endlink
    //! @author brownbp1, mendenjl
    //! @date Dec 11, 2018
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeMultiAlign :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! the matrix, which will hold the tanimoto coefficients that are found
      mutable double m_MultiAlignScore;

      //! filename for the first molecular ensemble
      util::ShPtr< command::ParameterInterface> m_InputMolSDF;

      util::ShPtr< command::ParameterInterface> m_InputPocketPDB;

      //! flag indicating the implementation of chemistry::ComformationComparisonInterface to use to compare conformers
      util::ShPtr< command::FlagInterface> m_ConformerComparerFlag;

      //! flag indication of first molecule to load from ensemble A
      util::ShPtr< command::FlagInterface> m_StartA;

      //! flag indication of last molecule to load from ensemble A
      util::ShPtr< command::FlagInterface> m_MaxMolsA;

      mutable chemistry::ConformationComparisonMultiAlign m_Comparers;
      mutable util::SiPtrVector< const chemistry::ConformationInterface> m_EnsembleA;

      //! ensemble sizes
      mutable size_t m_EnsembleASize;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      MoleculeMultiAlign();

    public:

      // instantiate enumerator for
      static const ApplicationType MoleculeMultiAlign_Instance;

      //! @brief Clone function
      //! @return pointer to new FoldProtein
       MoleculeMultiAlign *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

//    private:
//
//      //! @brief compare the molecules given by the indices in a vector
//      void RunThread( const size_t &THREAD_ID) const;

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

    }; //

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MOLECULE_MULTIALIGN_H_
