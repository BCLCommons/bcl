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

#ifndef BCL_RESTRAINT_ANALYZE_ATOM_DISTANCE_HEATMAP_H_
#define BCL_RESTRAINT_ANALYZE_ATOM_DISTANCE_HEATMAP_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_restraint_interface.h"
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AnalyzeAtomDistanceHeatmap
    //! @brief makes a heat map showing the frequency with which a distance occurs for each atom distance restraint
    //! @details can optionally also show in the heat map the restraint distance with upper and lower bounds
    //!
    //! @see @link example_restraint_analyze_atom_distance_heatmap.cpp @endlink
    //! @author alexanns
    //! @date Aug 1, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AnalyzeAtomDistanceHeatmap :
      public assemble::AnalyzeProteinEnsembleInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the postfix that will be appended to the filename that will hold this analysis
      std::string m_OutFilePostFix;

      //! the type of protein model data this atom distance restraint corresponds to
      util::Implementation< HandlerBase< util::ShPtrVector< AtomDistance> > > m_ProteinModelDataType;

      //! he minimal value representing the left boundary of the score histogram
      double m_HistogramMinimum;

      //! the width of one bin of the score histograms
      double m_HistogramBinSize;

      //! the number of bins in the score histograms
      size_t m_HistogramNumberOfBins;

      //! true if the atom distance restraint should be indicated on the heatmap
      bool m_ShowAtomDistance;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AnalyzeAtomDistanceHeatmap();

      //! @brief Clone function
      //! @return pointer to new AnalyzeAtomDistanceHeatmap
      AnalyzeAtomDistanceHeatmap *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives the string to append to the the end of a filename to identify this analysis
      //! @return gives the string to append to the the end of a filename to identify this analysis
      const std::string &GetOutFilePostfix() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
      //! @param ENSEMBLE the ensemble that will be analyzed
      //! @return string which has the analyzed information about the ensemble
      std::string operator()( const assemble::ProteinEnsemble &ENSEMBLE) const;

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief creates histograms of how frequently a distance is observed in models for each restraint
      //! @param DATA the restraints that will be used to get distances
      //! @param ENSEMBLE the ensemble of models that will be used to get distances
      //! @return vector of histograms - one histogram for each restraint
      storage::Vector< math::Histogram> GetDistanceHistograms
      (
        const util::ShPtrVector< AtomDistance> &DATA, const assemble::ProteinEnsemble &ENSEMBLE
      ) const;

      //! @brief creates the tics that will be used in the heatmap - converts each restraint into a string to use as tic
      //! @param DATA the restraints that will be used to create tics
      //! @return vector of strings which are the tics representing each restraint
      static storage::Vector< std::string> GetRestraintNameTics
      (
        const util::ShPtrVector< AtomDistance> &DATA
      );

      //! @brief creates gnuplot heat map object from the distance histograms and the restraint names as tics
      //! @param HISTOGRAMS the histograms that will be used to make the heat map
      //! @param TICS the names of the restraints
      //! @return heat map object which represents the distribution of distances for each restraint
      static math::GnuplotHeatmap GetHeatMap
      (
        const storage::Vector< math::Histogram> &HISTOGRAMS, const storage::Vector< std::string> &TICS
      );

      //! @brief gives the coordinates of boxes representing the atom distance restraint data
      //!        each atom distance restraint gives two boxes, one for the upper portion of the restraint, and
      //!        one for the lower portion of the restraint. The actual distance is given by the junction of the two
      //!        boxes
      //! @param DATA the restraints that will be used to create the boxes
      //! @return vector of [(x,y),(x,y)] coordinates - the lower left corner and upper right corner of the box,
      //!         respectively
      static storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > >
      GetDistanceBoxes( const util::ShPtrVector< AtomDistance> &DATA);

    }; // class AnalyzeAtomDistanceHeatmap

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_ANALYZE_ATOM_DISTANCE_HEATMAP_H_
