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

#ifndef BCL_ASSEMBLE_PROTEIN_ENSEMBLE_H_
#define BCL_ASSEMBLE_PROTEIN_ENSEMBLE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "find/bcl_find.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "restraint/bcl_restraint.fwd.hh"
#include "sspred/bcl_sspred.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_collector_common_aa.h"
#include "bcl_assemble_domain.h"
#include "bcl_assemble_protein_model.h"
#include "bcl_assemble_sse_geometry.h"
#include "biol/bcl_biol_atom.h"
#include "math/bcl_math_limits.h"
#include "math/bcl_math_running_average_sd.h"
#include "util/bcl_util_sh_ptr_list.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinEnsemble
    //! @brief represents a collection (ensemble) of proteins
    //! @details provides basic functionality useful for an ensemble of proteins
    //!
    //! @see @link example_assemble_protein_ensemble.cpp @endlink
    //! @author alexanns
    //! @date Feb 11, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinEnsemble :
      public ProteinModel
    {

    private:

    //////////
    // data //
    //////////

      //! the protein models comprising this ensemble
      util::ShPtrVector< ProteinModel> m_Ensemble;

    public:

      //! typedef for iterator
      typedef std::vector< util::ShPtr< ProteinModel> >::iterator       iterator;

      //! typedef for const_iterator
      typedef std::vector< util::ShPtr< ProteinModel> >::const_iterator const_iterator;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! command line flag for reading list of protein models
      static util::ShPtr< command::FlagInterface> &GetFlagPDBList();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ProteinEnsemble();

      //! @brief constructor
      //! @param MODEL model the ensemble will be created from
      ProteinEnsemble( const util::ShPtr< ProteinModel> &MODEL);

      //! @brief constructor
      //! @param MODEL model the ensemble will be created from
      ProteinEnsemble( const ProteinModel &MODEL);

      //! @brief constructor from range of iterators
      //! @param FIRST iterator to first element to be copied
      //! @param LAST iterator to one past the last element to be copied
      template< typename t_InputIterator>
      ProteinEnsemble( const t_InputIterator &FIRST, const t_InputIterator &LAST);

      //! @brief constructor taking filename with list of pdbs in it
      //! @param FILENAME file with the list of pdbs in it
      //! @param COLUMN the column in the file that has the pdbs in it
      //! @param AA_CLASS aa class used for protein models
      //! @param PREFIX prefix to add to all filenames in the file
      ProteinEnsemble
      (
        const std::string &FILENAME,
        const size_t COLUMN,
        const biol::AAClass &AA_CLASS,
        const std::string &PREFIX = std::string(),
        const bool &STATUS_MESSAGE = false,
        const size_t &INPUT_START = 0,
        const size_t &INPUT_MAX = math::GetHighestBoundedValue< size_t>()
      );

      //! @brief Clone function
      //! @return pointer to new MutateDistanceDistribution
      ProteinEnsemble *Clone() const;

      //! @brief hard copy constructor
      //! @return a ProteinModel with chains hard copied from that model
      virtual ProteinEnsemble *HardCopy() const;

      //! @brief empty copy constructor
      //! @return a ProteinEnsemble that is empty
      virtual ProteinEnsemble *Empty() const;

      //! @brief destructor
      virtual ~ProteinEnsemble();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      util::ShPtrVector< ProteinModel> &GetEnsembleData()
      {
        return m_Ensemble;
      }

      const util::ShPtrVector< ProteinModel> &GetEnsembleData() const
      {
        return m_Ensemble;
      }

      //! @brief returns size of the container
      //! @return size, i.e. number of elements stored
      size_t GetSize() const;

      //! @brief return iterator on begin
      //! @return iterator pointing to the beginning of the container, i.e. the first element
      iterator Begin();

      //! @brief return const_iterator on begin
      //! @return const_iterator pointing to the beginning of the container, i.e. the first element
      const_iterator Begin() const;

      //! @brief return iterator on end
      //! @return iterator pointing to the end of the container, i.e. behind the last element
      iterator End();

      //! @brief return const_iterator on end
      //! @return const_iterator pointing to the end of the container, i.e. behind the last element
      const_iterator End() const;

      //! @brief sets the sse pool within the protein model data
      //! @param SSE_POOL the sse pool that will be set in the protein model data
      virtual void SetSSEPoolData( const util::ShPtr< SSEPool> &SSE_POOL);

      //! @brief set identifiers of each model in the ensemble
      //! @param IDENTIFIERS the identifiers in the order of the models that they are assigned
      //! @return true if setting was successful
      bool SetIdentifiers( const storage::Vector< std::string> &IDENTIFIERS);

    ////////////////
    // operations //
    ////////////////

      //! @brief checks whether container is empty
      //! @return if the container is empty or not
      bool IsEmpty() const;

      //! @brief insert ELEMENT into the container
      //! @param ELEMENT an object of t_DataType that is inserted
      void InsertElement( const util::ShPtr< ProteinModel> &ELEMENT);

      //! @brief insert ELEMENT into the container
      //! @param POS the position where the element should be inserted
      //! @param ELEMENT an object of t_DataType that is inserted
      void InsertElement( const size_t POS, const util::ShPtr< ProteinModel> &ELEMENT);

      //! @brief delete single element at ITR
      //! @param ITR iterator pointing to the element that will be destroyed
      void RemoveElement( iterator ITR);

      //! @brief orders the models in the ensemble according to a comparison function
      //! @param COMP comparison function object returning true if first argument comes before the second argument
      template< typename t_Comparison>
      void Sort( const t_Comparison &COMP);

      //! @brief determines which proteins in this ensemble are not contained in the other ensemble
      //!        both ensembles must be sorted!
      //! @param OTHER_ENSEMBLE the protein ensemble that will be subtracted from this ensemble
      //! @param COMP comparison function object to determine equivalence
      //! @return ProteinEnsemble which has the proteins contained in this ensemble that are not in OTHER_ENSEMBLE
      template< typename t_Comparison>
      ProteinEnsemble Difference( ProteinEnsemble &OTHER_ENSEMBLE, const t_Comparison &COMP) const;

      //! @brief gives the mean and std dev of distances for a data pair for this ensemble
      //! @param DATA_PAIR the data pair that indicates the distance that should be calculated
      //! @return the mean and std dev of all pairwise distances
      math::RunningAverageSD< double> GetDistanceStatistics( const restraint::DataPairwise &DATA_PAIR) const;

      //! @brief gives vector of coordinates located from each of the models in the ensemble
      //! @param COORD_LOCATOR method for locating coordinates in each model
      //! @return vector of coordinates located from each of the models in the ensemble
      storage::Vector< linal::Vector3D> GetCoordinates
      (
        const find::LocatorInterface< linal::Vector3D, ProteinModel> &COORD_LOCATOR
      ) const;

      //! @brief gives a vector of distances calculated within each of the models for a given data pair
      //! @param DATA_PAIR the data pair that indicates the distance that should be calculated
      //! @return vector of doubles which are the distances calculated in each of the models of the ensemble
      storage::Vector< double> GetDistances( const restraint::DataPairwise &DATA_PAIR) const;

      //! @brief gives a vector of distance changes for a data pair for this ensemble versus a given ensemble.
      //!        all pairwise distance changes are calculated and provided
      //! @param DATA_PAIR the data pair that indicates the distance that should be calculated
      //! @param OTHER_ENSEMBLE the ensemble which provides the second state to get distance changes from
      //! @return vector of doubles which are all pairwise distance changes going from this ensemble to OTHER_ENSEMBLE
      //!         for the given data pair
      storage::Vector< double> GetDistanceChanges
      (
        const restraint::DataPairwise &DATA_PAIR,
        const ProteinEnsemble &OTHER_ENSEMBLE
      ) const;

      //! @brief gives the mean and std dev of distance changes for a data pair for this ensemble versus a given
      //!        ensemble. all pairwise distance changes are calculated and provided
      //! @param DATA_PAIR the data pair that indicates the distance that should be calculated
      //! @param OTHER_ENSEMBLE the ensemble which provides the second state to get distance changes from
      //! @return the mean and std dev of all pairwise distance changes going from this ensemble to OTHER_ENSEMBLE
      math::RunningAverageSD< double> GetDistanceChangesMeanSD
      (
        const restraint::DataPairwise &DATA_PAIR,
        const ProteinEnsemble &OTHER_ENSEMBLE
      ) const;

      //! @brief gives all mean and std dev of distance changes for a data set for this ensemble versus a given
      //!        ensemble. all pairwise distance changes are calculated and provided
      //! @param DATA_SET the data set that indicates the distances that should be calculated
      //! @param OTHER_ENSEMBLE the ensemble which provides the second state to get distance changes from
      //! @return the mean and std dev of all pairwise distance changes going from this ensemble to OTHER_ENSEMBLE
      std::multimap< math::RunningAverageSD< double>, restraint::DataPairwise, math::LessThanAbsoluteMean>
      GetDistanceChangesMeanSD
      (
        const restraint::DataSetPairwise &DATA_SET, const ProteinEnsemble &OTHER_ENSEMBLE
      ) const;

      //! @brief removes all models from the ensemble
      void Reset();

      //! @brief read given ssmethods predictions
      //! @param SS_METHODS the methods to read for each protein in the ensemble
      //! @return true if reading was successful
      bool ReadSSPredictions( const storage::Set< sspred::Method> &SS_METHODS);

      //! @brief gives all of the pdb file names of the proteins in the ensemble
      //! @return vector of strings which are the pdb names of the proteins in the ensemble
      storage::Vector< std::string> GetPDBNames() const;

      //! @brief gives a format object which can be used to format strings based on the longest pdb name
      //! @return format object which can be used to format strings to the same length as the longest pdb name
      util::Format GetNameFormatter() const;

      //! @brief gives the conformations that are in this protein
      //! @return protein ensemble which is the conformations making up this protein
      virtual const ProteinEnsemble &GetConformationalEnsemble() const;

      //! @brief sets the conformations that are in this protein
      //! @param ENSEMBLE protein ensemble which is the conformations making up this protein
      void SetConformationalEnsemble( const ProteinEnsemble &ENSEMBLE);

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

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief opens a file and reads a list of pdbs and creates an ensemble out of them
      //! @param FILENAME file with the list of pdbs in it
      //! @param COLUMN the column in the file that has the pdbs in it - start at 0
      //! @param PREFIX prefix to add to all filenames in the file
      static ProteinEnsemble GetEnsembleFromFile
      (
        const std::string &FILENAME,
        const size_t COLUMN,
        const biol::AAClass &AA_CLASS,
        const std::string &PREFIX = std::string(),
        const bool &STATUS_MESSAGE = false,
        const size_t &INPUT_START = 0,
        const size_t &INPUT_MAX = math::GetHighestBoundedValue< size_t>()
      );

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment operator
      //! @param PROTEIN_MODEL ProteinModel to be copied
      //! @return This model after all members are assigned to values from PROTEIN_MODEL
      ProteinEnsemble &operator =( const ProteinModel &PROTEIN_MODEL)
      {
        SetConformationalEnsemble( PROTEIN_MODEL.GetConformationalEnsemble());

        ProteinModel::operator =( PROTEIN_MODEL);

        return *this;
      }

    }; // class ProteinEnsemble

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from range of iterators
    //! @param FIRST iterator to first element to be copied
    //! @param LAST iterator to one past the last element to be copied
    template< typename t_InputIterator>
    ProteinEnsemble::ProteinEnsemble( const t_InputIterator &FIRST, const t_InputIterator &LAST) :
      m_Ensemble( FIRST, LAST)
    {
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief orders the models in the ensemble according to a comparison function
    //! @param COMP comparison function object returning true if first argument comes before the second argument
    template< typename t_Comparison> void ProteinEnsemble::Sort( const t_Comparison &COMP)
    {
      m_Ensemble.Sort( COMP);
    }

    //! @brief determines which proteins in this ensemble are not contained in the other ensemble
    //!        this ensemble must be sorted before calling this!
    //! @param OTHER_ENSEMBLE the protein ensemble that will be subtracted from this ensemble
    //! @param COMP comparison function object to determine equivalence
    //! @return ProteinEnsemble which has the proteins contained in this ensemble that are not in OTHER_ENSEMBLE
    template< typename t_Comparison> ProteinEnsemble ProteinEnsemble::Difference
    (
      ProteinEnsemble &OTHER_ENSEMBLE, const t_Comparison &COMP
    ) const
    {
      // sort "OTHER_ENSEMBLE"
      OTHER_ENSEMBLE.Sort( COMP);

      // ProteinEnsemble that will hold the models that are in this ensemble but not in "OTHER_ENSEMBLE"
      ProteinEnsemble difference_models;

      // find the models that are in this ensemble but not in "OTHER_ENSEMBLE"
      // put the models in "difference_models"
      std::set_difference
      (
        Begin(), End(), OTHER_ENSEMBLE.Begin(), OTHER_ENSEMBLE.End(),
        std::inserter( difference_models.m_Ensemble.InternalData(), difference_models.Begin()), COMP
      );

      // return the new ensemble holding the models that are in this ensemble but not in "OTHER_ENSEMBLE"
      return difference_models;
    }

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_PROTEIN_ENSEMBLE_H_
