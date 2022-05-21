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

#ifndef BCL_RESTRAINT_DATA_PAIRWISE_H_
#define BCL_RESTRAINT_DATA_PAIRWISE_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom_coordinates_interface.h"
#include "math/bcl_math_running_min_max.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataPairwise
    //! @brief indicates a piece of data that needs two points
    //! @details Two data points that are needed to indicate a meaningful piece of data. The data points are sorted
    //!          always by sequence.
    //!
    //! @see @link example_restraint_data_pairwise.cpp @endlink
    //! @author alexanns
    //! @date May 6, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DataPairwise :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the first point in the data pair
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface> m_First;

      //! the second point in the data pair
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface> m_Second;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DataPairwise();

      //! @brief constructor taking members
      //! @param LOCATOR_A first locator
      //! @param LOCATOR_B second locator
      DataPairwise
      (
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_A,
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_B
      );

      //! @brief Clone function
      //! @return pointer to new DataPairwise
      DataPairwise *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives formatted string describing the data pair
      //! @return formatted string describing the data pair
      std::string GetIdentification() const;

      //! @brief reads formatted string describing the locator
      //! @return formatted string describing the locator
      std::istream &ReadIdentification( std::istream &ISTREAM);

      //! @brief gives the first point in the data pair
      //! @return the first point in the data pair
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &First() const;

      //! @brief gives the second point in the data pair
      //! @return the second point in the data pair
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &Second() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief sets the two data points involved in this data pair
      //! @param ATOM_A the  first data point that will be set
      //! @param ATOM_B the second data point that will be set
      void Set
      (
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &ATOM_A,
        const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &ATOM_B
      );

      //! @brief indicates whether or not the two data points have been set
      //! @return boolean true if the data pair has been successfully set - false otherwise
      bool IsSet() const;

      //! @brief calculates the euclidian distance indicated by a protein model
      //! @param MODEL the model from which the distance will be calculated
      //! @return double which is the distance between the two coordinates indicated by this data pair
      double EuclidianDistance( const assemble::ProteinModel &MODEL) const;

      //! @brief calculates statistics for the distance indicated by this across an ensemble of models
      //! @param ENSEMBLE the ensemble over which distances and statistics will be calculated
      //! @return pair of RunningAverageSD< double> and RunningMinMax< double> indicated distance statistics over ENSEMBLE
      storage::Pair< math::RunningAverageSD< double>, math::RunningMinMax< double> >
      EuclidianDistance( const assemble::ProteinEnsemble &ENSEMBLE) const;

      //! @brief calculates the sequence separation between the two data points in this data pair
      //! @return size_t which is the sequence separation between the two data points in this data pair
      size_t SequenceSeparation() const;

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
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief reads formatted string describing the data pair
      //! @return DataPairwise read from formatted string describing the data pair
      DataPairwise ReadIdentification( std::istream &ISTREAM) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief calculates sequence separation between two residues indicated by chain and seq ids
      //! @param CHAIN_ID_A the chain id of the first residue
      //! @param SEQ_ID_A the seq id of the first residue
      //! @param CHAIN_ID_B the chain id of the second residue
      //! @param SEQ_ID_B the seq id of the second residue
      static size_t CalculateSequenceSeparation
      (
        const char CHAIN_ID_A, const int SEQ_ID_A, const char CHAIN_ID_B, const int SEQ_ID_B
      );

    }; // class DataPairwise

    //! @brief less than operator for comparing two DataPairwise
    //! @param LHS the first DataPairwise which will be compared against the second DataPairwise
    //! @param RHS the second DataPairwise which will be compared against the first DataPairwise
    //! @return boolean true if LHS is less than RHS - false otherwise
    bool operator <( const DataPairwise &LHS, const DataPairwise &RHS);

    //! @brief writes pymol script formatted data to a stream that can show distances in pymol
    //! @param OSTREAM the stream which will write the data
    //! @param DISTANCES the data that will be written along with their distance and the desired color of their line
    //! @return ostream the stream that wrote the data distance information
    std::ostream &ShowDistancesInPymol
    (
      std::ostream &OSTREAM,
      const storage::List< storage::Triplet< DataPairwise, double, linal::Vector3D> > &DISTANCES
    );

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_DATA_PAIRWISE_H_
