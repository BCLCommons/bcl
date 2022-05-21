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

#ifndef BCL_APP_MOLECULE_COORDINATES_H_
#define BCL_APP_MOLECULE_COORDINATES_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_types.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_interface.h"
#include "math/bcl_math_running_average_sd.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_vector_nd.h"

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeCoordinates
    //! @brief Application for processing molecular coordinates.  Current options include analysis and histograms of
    //!         bonds and angles in small molecule ensembles, also recentering molecules
    //!
    //! @see @link example_app_molecule_coordinates.cpp @endlink
    //! @author mendenjl, brownbp1
    //! @date Mar 7, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeCoordinates :
      public InterfaceRelease
    {
    private:

    //////////
    // data //
    //////////

      //! Dihedral bin size (in degrees)
      util::ShPtr< command::FlagInterface> m_DihedralBinSizeInDegrees;

      //! Bond angle bin size (in degrees)
      util::ShPtr< command::FlagInterface> m_BondAngleBinSizeInDegrees;

      //! Whether to compute statistics
      util::ShPtr< command::FlagInterface> m_StatisticsFlag;

      //! Whether to compute dihedral scores
      util::ShPtr< command::FlagInterface> m_DihedralScoresFlag;

      //! Whether to compute amide deviations (in degrees)
      util::ShPtr< command::FlagInterface> m_AmideBondDeviationsFlag;

      //! Whether to compute amide penalties given standard tolerance
      //! values of 10, 15, and 25 degrees with corresponding penalties
      //! of 0.01, 0.1, and 1.0.
      util::ShPtr< command::FlagInterface> m_AmideBondPenaltiesFlag;

      //! Whether to compute per-atom-pair clash scores
      util::ShPtr< command::FlagInterface> m_ClashScoresFlag;

      // whether to recenter the molecules
      util::ShPtr< command::FlagInterface> m_RecenterMoleculesFlag;

      // whether to output centroid coordinates
      util::ShPtr< command::FlagInterface> m_MoleculeCentroidFlag;

      //! Output filename base
      util::ShPtr< command::FlagInterface> m_OutputFilenameBase;

      //! dihedral angle histogram, keyed by a string containing the atom types and bond types involved
      mutable storage::Map< std::string, math::Histogram> m_DihedralHistogram;

      //! bond angle histogram, keyed by a string containing the atom types and bond types involved
      mutable storage::Map< std::string, math::Histogram> m_BondAngleHistogram;

      //! bond angle histogram, keyed by a string containing the atom types and bond types involved
      mutable storage::Map< std::string, math::RunningAverageSD< double> > m_BondAngleStats;

      //! bond length histogram, keyed by a string containing the atom types and bond types involved
      mutable storage::Map< std::string, math::RunningAverageSD< double> > m_BondLengths;

      //! bond length histogram, keyed by a string containing the atom types and bond types involved
      mutable storage::Map< std::string, math::RunningAverageSD< double> > m_BondDihedrals;

      //! dihedral angle scores
//      mutable storage::Vector<>

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      MoleculeCoordinates();

    public:

      // instantiate enumerator for MoleculeCoordinates class
      static const ApplicationType MoleculeCoordinates_Instance;

      //! @brief Clone function
      //! @return pointer to new MoleculeCoordinates
      MoleculeCoordinates *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const;

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

      //! @brief write histogram given in a map with key identifier to a stream
      //! @param OUTPUT output stream to write to
      //! @param HISTOGRAMS histogram in map with key identifier
      //! @param REFLECTED true if the degrees are reflective (e.g. dihedral bins)
      static void WriteHistogram
      (
        std::ostream &OUTPUT,
        const storage::Map< std::string, math::Histogram> &HISTOGRAMS,
        const double &BIN_SIZE,
        const bool &REFLECTED = false
      );

      //! @brief write statistics given in a map with key identifier to a stream
      //! @param OUTPUT output stream to write to
      //! @param STATS statistics in map with key identifier
      //! @param BIN_SIZE bin size for the statistics
      static void WriteStatistics
      (
        std::ostream &OUTPUT,
        const storage::Map< std::string, math::RunningAverageSD< double> > &STATS
      );

      //! @brief write centroid coordinates to a stream
      //! @param OUTPUT output stream to write to
      //! @param INDEX molecule feed index
      //! @param CENTROID coordinate vector for molecule center
      static void WriteCentroid
      (
        std::ostream &OUTPUT,
        size_t INDEX,
        const linal::Vector3D CENTROID
      );

      //! @brief write dihedral scores to output
      //! @param OUTPUT output stream to write to
      //! @param INDEX molecule feed index
      //! @param SCORES scores and corresponding atom indices
      static void WriteScores
      (
        std::ostream &OUTPUT,
        size_t INDEX,
        const storage::Vector< storage::Triplet< size_t, size_t, double> > SCORES
      );

      //! @brief returns a simple name for a bond type given as a size_t into a string
      //! @returns one of XXXXXX, Aromatic, [1-3]x(Chain|Ring) depending on the bond type
      static std::string ConvertBondTypeIDToString( const size_t &BOND_TYPE_ID);

      //! @brief returns a simple name for a bond type given as a size_t into a string
      //! @returns one of XXXXXX, Aromatic, Single, Double, or Triple depending on the bond type
      static std::string ConvertBondTypeIDToStringIgnoreTopology( const size_t &BOND_TYPE_ID);

      //! @brief add the data from a small molecule to the histograms
      //! @param MOLECULE the small molecule whose data to add
      void AddDataToHistograms( const chemistry::ConformationInterface &MOLECULE) const;

      //! @brief compute the dihedral angle scores for the molecule
      //! @param MOLECULE the small molecule whose dihedral scores will be computed
      //! @return a vector of scores and the atom indices of the central dihedral bond atoms
      static storage::Vector< storage::Triplet< size_t, size_t, double> > ComputeDihedralScores
      (
        const chemistry::ConformationInterface &MOLECULE
      );

    public:
      //! @brief check whether the given vector nd refers to a non-gasteiger atom type
      //! @param TYPE vector containing alternating atom/bond types, starting with an atom type
      //! @return true if all the atom types are proper gasteiger types
      template< unsigned int t_N>
      static bool ContainsOnlyGasteigerAtomTypes( const storage::VectorND< t_N, size_t> &TYPE)
      {
        // determine the last atom type that is a gasteiger type
        static size_t n_gasteiger_types( chemistry::AtomTypes::GetNumberGasteigerTypes());

        for( unsigned int i( 0); i < t_N; i += 2)
        {
          if( TYPE( i) >= n_gasteiger_types)
          {
            return false;
          }
        }
        return true;
      }

      //! @brief check whether the given vector nd refers to a non-gasteiger atom type
      //! @param TYPE vector containing alternating atom/bond types, starting with an atom type
      //! @return true if all the atom types are proper gasteiger types
      template< unsigned int t_N>
      static std::string ConvertAtomBondTypeIntoString( const storage::VectorND< t_N, size_t> &TYPE)
      {
        std::ostringstream oss;
        oss << chemistry::AtomType( TYPE( 0)).GetName() << ' ';
        for( unsigned int i( 1); i < t_N; i += 2)
        {
          oss << ConvertBondTypeIDToString( TYPE( i)) << ' ';
          oss << chemistry::AtomType( TYPE( i + 1)).GetName() << ' ';
        }
        return oss.str();
      }

    }; // MoleculeCoordinates

  } // namespace app
} // namespace bcl

#endif // BCL_APP_MOLECULE_COORDINATES_H_
