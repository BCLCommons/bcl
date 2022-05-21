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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#ifndef BCL_APP_CONFORMER_GENERATOR_H_
#define BCL_APP_CONFORMER_GENERATOR_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"
#include "descriptor/bcl_descriptor.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_comparison_interface.h"
#include "io/bcl_io_directory.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConformerGenerator
    //! @brief Application that samples molecule conformations
    //!
    //! @see @link example_app_conformer_generator.cpp @endlink
    //! @author kothiwsk
    //! @date Aug 29, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConformerGenerator :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! file containing the fragment rotamer library
      util::ShPtr< command::FlagInterface>                                        m_RotamerLibraryFlag;

      //! file containing native conformation that generated conformers will be compared against
      util::ShPtr< command::FlagInterface>                                        m_NativeEnsembleFileFlag;

      //! flag indicating the implementation of chemistry::ComformationComparisonInterface to use to compare conformers
      util::ShPtr< command::FlagInterface>                                        m_ConformerComparerFlag;

      //! conformers will only be considered a novel if conformation comparison value is greater that
      util::ShPtr< command::FlagInterface>                                        m_MaxIterations;

      //! prefix for output filename that stores conformations of a particular molecule
      util::ShPtr< command::FlagInterface>                                        m_OutputLigandData;

      //! flag to specify number of top models that are required
      util::ShPtr< command::FlagInterface>                                        m_TopModels;

      //! prefix for output filename for molecules whose 3D structures could not be generated
      util::ShPtr< command::FlagInterface>                                        m_OutputFailed3D;

      //! prefix for file containing rmsd and score data
      util::ShPtr< command::FlagInterface>                                        m_RmsdScorePrefix;

      //! prefix for file containing rmsd and score data
      util::ShPtr< command::FlagInterface>                                        m_ConformersSingleFile;

      util::ShPtr< command::FlagInterface>                                        m_ConformersSeparateFile;

      util::ShPtr< command::FlagInterface>                                        m_FixedNumberConformers;

      util::ShPtr< command::FlagInterface>                                        m_ChangeChirality;

      util::ShPtr< command::FlagInterface>                                        m_Generate3DCoordinates;

      util::ShPtr< command::FlagInterface>                                        m_RandomSelection;

      util::ShPtr< command::FlagInterface>                                        m_Cluster;      

      util::ShPtr< command::FlagInterface>                                        m_NoCluster;

      util::ShPtr< command::FlagInterface>                                        m_RandomDihedralMutateWeight;

      util::ShPtr< command::FlagInterface>                                        m_SkipDihedralSamplingFlag;

      util::ShPtr< command::FlagInterface>                                        m_SkipBondAnglesFlag;

      util::ShPtr< command::FlagInterface>                                        m_SkipRingSamplingFlag;

      util::ShPtr< command::FlagInterface>                                        m_MaxClashResolutionIterationsFlag;

      util::ShPtr< command::FlagInterface>                                        m_MaxClashToleranceFlag;

      //! obtains a dissimilarity
      mutable util::Implementation< chemistry::ConformationComparisonInterface>   m_Comparer;

      //! class that searches fragments for a molecule of interest from fragment library
      mutable util::ShPtr< chemistry::SampleConformations>                        m_SampleConformations;

      //! class that maps fragments to a molecule of interest
      mutable util::ShPtr< chemistry::SmallMoleculeFragmentMapping>               m_FragmentMapping;

      //! container to store native conformations
      mutable storage::Vector< chemistry::FragmentComplete>                       m_NativeEnsemble;

      //! the output file, which writes out the conformers as they are generated
      mutable io::OFStream                                                        m_OutputSatisfied;

      //! output file which writes out rmsd
      mutable io::OFStream                                                        m_RmsdFile;

      //! output file which writes out rmsd
      mutable io::OFStream                                                        m_DihedralBinFile;

      //! output file which writes out score
      mutable io::OFStream                                                        m_ScoreFile;

      //! output file which writes out score
      mutable io::OFStream                                                        m_RmsdScore;

      //! directory where conformations are written
      mutable io::Directory                                                       m_Directory;

      //! the output file, which writes out the conformers as they are generated
      mutable io::OFStream                                                        m_OutputConformers;

      //! the output file, which writes out the conformers as they are generated
      mutable io::OFStream                                                        m_Failed3DStream;

      //! bin size to be used for comparing conformations
      mutable double                                                              m_ConformationComparison;

      //! bin size to be used for comparing conformations
      mutable double                                                              m_Diversity;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      ConformerGenerator();

      //! copy constructor
      ConformerGenerator( const ConformerGenerator &APP);

    public:

      // instantiate enumerator for PrepareSmallMoleculeEnsemble class
      static const ApplicationType ConformerGenerator_Instance;

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      ConformerGenerator *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the full description
      //! @return the full description
      std::string GetDescription() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetWebText() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

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

      //! @brief function that writes rmsd and score data to require file stream
      //! @param MOLECULE molecule of interest whose conformations need to be sampled
      //! @param MOLECULE_INDEX index of molecule of interest in ensemble
      //! @param ENSEMBLE ensemble that contains conformations of molecule of interest
      void RmsdScoreFile
      (
        const chemistry::FragmentComplete &MOLECULE,
        const size_t &MOLECULE_INDEX,
        const chemistry::FragmentEnsemble &ENSEMBLE
      ) const;

      //! @brief function that creates pymol representation of fragments and conformations of molecule of interest
      //! @param MOLECULE molecule of interest whose conformations need to be sampled
      //! @param MOLECULE_INDEX index of molecule of interest in ensemble
      //! @param ENSEMBLE ensemble that contains conformations of molecule of interest
      //! @param FRAGMENT_ISO fragments that are part of molecule of interest
      void PymolRepresentation
      (
        const chemistry::FragmentComplete &MOLECULE,
        const size_t &MOLECULE_INDEX,
        const chemistry::FragmentEnsemble &ENSEMBLE,
        const chemistry::FragmentEnsemble &FRAGMENT_ENSEMBLE
      ) const;

      //! @brief controls output from the app of interest
      //! @param MOLECULES the molecule of interest whose conformations have to be sampled
      //! @param ENSEMBLE ensemble of conformations that were sampled for the molecule of interest
      //! @param NUMBER_OF_CONFORMATIONS number of conformations desired
      //! @param MOLECULE_INDEX index of the molecule in the ensemble
      void Output
      (
        const chemistry::FragmentComplete MOLECULE,
        const chemistry::FragmentEnsemble &ENSEMBLE,
        const size_t NUMBER_OF_CONFORMATIONS,
        const size_t MOLECULE_INDEX
      ) const;

      //! @brief generates conformations for a molecule of interest
      //! @param MOLECULES the molecule of interest whose conformations have to be sampled
      //! @param MOLECULE_INDEX index of the molecule in the ensemble
      void ConformationSampling
      (
        const chemistry::FragmentComplete &MOLECULE,
        const size_t &MOLECULE_INDEX
      ) const;

    }; // ConformerGeneration

  } // namespace app
} // namespace bcl
#endif // BCL_APP_CONFORMER_GENERATOR_H_
