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

#ifndef BCL_SCORESTAT_STRAND_ALIGNMENT_H_
#define BCL_SCORESTAT_STRAND_ALIGNMENT_H_

// include the namespace header
#include "bcl_scorestat.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "math/bcl_math_histogram.h"
#include "storage/bcl_storage_table.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
    /////////////////////////////////////////////////////////////////////
    //!
    //! @class StrandAlignment
    //! @brief extracts strand alignment statistics from protein models
    //!
    //! @see @link example_strand_alignment_statistics.cpp @endlink
    //! @author lib14
    //! @date Dec 1, 2014
    //////////////////////////////////////////////////////////////////////

    class BCL_API StrandAlignment :
      public assemble::AnalyzeProteinEnsembleInterface
    {

    public:
        //! output options
        enum OutputOption
        {
          e_Table,
          s_NumberOutputOptions
        };

        //! @brief OutputOption as string
        //! @param OUTPUT_OPTION the OutputOption
        //! @return the string for the OutputOption
        static const std::string &GetOutputOptionName( const OutputOption &OUTPUT_OPTION);

        //! @brief Output filename as string
        //! @param OutputOption the desired Output Type
        //! @return the string for the output file extension
        static const std::string &GetOutputFileName( const OutputOption &OUTPUT_OPTION);

        //! @brief OutputOptionEnum enum I/O helper
        typedef util::WrapperEnum< OutputOption, &GetOutputOptionName, s_NumberOutputOptions> OutputOptionEnum;

    private:

    //////////
    // data //
    //////////

        //! output options
        OutputOptionEnum m_OutputOption;

        //! cutoff for the distance between Carbonyl-Oxygen and Amide-Nitrogen
        double m_HydrogenBondOHCutoff;

        //! path where to write the pdbs to visualize measurements
        std::string m_VisualizationOutputPath;

        //! chain ids
        std::string m_ChainIds;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      StrandAlignment();

      //! @brief virtual copy constructor
      StrandAlignment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as constant reference to std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives the string to append to the the end of a filename to identify this analysis
      //! @return gives the string to append to the the end of a filename to identify this analysis
      const std::string &GetOutFilePostfix() const;

      //! @brief gets the cutoff for distance between Carbonyl-Oxygen and Amide_Nitrogen
      //! @return the cutoff for distance between Carbonyl-Oxygen and Amide_Nitrogen
      const double &GetHydrogenBondOHCutoff() const;

      //! @brief gets the path where to write the pdbs to visualize measurements
      //! @return the path where to write the pdbs to visualize measurements
      const std::string &GetVisualizationOutputPath() const;

      //! @brief gets chain ids
      //! @return chain ids
      const std::string &GetChainIds() const;

      //! @brief gets the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
      //! @param ENSEMBLE the protein ensemble that will be analyzed
      //! @return string which has the analyzed information about the ensemble
      std::string operator()( const assemble::ProteinEnsemble &ENSEMBLE) const;

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief computes angle between the given three atoms; not using the linal functions b/c of Normalize() call
      //! @param ATOM_A the atom for calculating vector B->A
      //! @param ATOM_B the vertex
      //! @param ATOM_C the atom for calculating vector B->C
      //! @return the angle in radians
      double ComputeAngle( const biol::Atom &ATOM_A, const biol::Atom &ATOM_B, const biol::Atom &ATOM_C) const;

      //! @brief computes the H-bond distance between C-O and H-N, returns undefined if H coordinates are not defined
      //! @brief this function use an oriented H-bond definition: C-O ---> H-N
      //! @param AA_A the first AA ShPtr, the O is taken from this
      //! @param AA_B the second AA ShPtr, the H is taken from this
      //! @return the distance; undefined if a shptr is undefined or if no coordinates for the H are present
      double ComputeHBondDistance( const util::ShPtr< biol::AABase> &AA_A, const util::ShPtr< biol::AABase> &AA_B) const;

      //! @brief computes a list of values for a single H-bond between the two given AAs
      //! @brief H-bond are defined as oriented: C-O ---> H-N, so the C-O is from AA_A, the H-N is from AA_B
      //! @param AA_A the first AA ShPtr, the C-O is taken from this
      //! @param AA_B the second AA ShPtr, the H-N is taken from this
      //! @param ORIENTATION orientation of the sheet
      //! @return a list of vectors of H-bond specific values identical to the final table columns
      storage::List< storage::Vector< double> > ComputeStrandData
      (
        const util::ShPtr< biol::AABase> &AA_A,
        const util::ShPtr< biol::AABase> &AA_B,
        const assemble::SSEGeometryPacking::Orientation &ORIENTATION
      ) const;

      //! @brief collect all atoms that are part of an H-bond; use an oriented H-bond definition: C-O ---> H-N
      //! @param AA_A the first AA ShPtr, the C-O is taken from this
      //! @param AA_B the second AA ShPtr, the H-N is taken from this
      //! @return the list of atoms of this H-bond
      storage::List< biol::Atom> CollectStrandAtoms
      (
        const util::ShPtr< biol::AABase> &AA_A, const util::ShPtr< biol::AABase> &AA_B
      ) const;

      //! @brief selects the shortest possible H-bond (below a cutoff) for each AA from a given set of AA pairs
      //! @param STRAND_A first strand, one AA in each AA pair comes from this strand
      //! @param STRAND_B second strand, one AA in each AA pair comes from this strand
      //! @param POSSIBLE_HBONDS map of possible H-bonds (pairs of AA seqids)
      //! @return a list of pairs of AA Shptr that are the shortest H-bonds below a cutoff
      storage::List< storage::Pair< util::ShPtr< biol::AABase>, util::ShPtr< biol::AABase> > > SelectShortestHBonds
      (
        const assemble::SSE &STRAND_A,
        const assemble::SSE &STRAND_B,
        const std::multimap< int, int> &POSSIBLE_HBONDS
      ) const;

    }; // end class StrandAlignment
  } // namespace scorestat
} // namespace bcl

#endif // BCL_SCORESTAT_STRAND_ALIGNMENT_H_
