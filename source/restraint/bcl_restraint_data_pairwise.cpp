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
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "restraint/bcl_restraint_data_pairwise.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_string_replacement.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DataPairwise::s_Instance
    (
      GetObjectInstances().AddInstance( new DataPairwise())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DataPairwise::DataPairwise() :
      m_First(),
      m_Second()
    {
    }

    //! @brief constructor taking members
    //! @param LOCATOR_A first locator
    //! @param LOCATOR_B second locator
    DataPairwise::DataPairwise
    (
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_A,
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &LOCATOR_B
    ) :
      m_First(),
      m_Second()
    {
      Set( LOCATOR_A, LOCATOR_B);
    }

    //! @brief Clone function
    //! @return pointer to new DataPairwise
    DataPairwise *DataPairwise::Clone() const
    {
      return new DataPairwise( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataPairwise::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives formatted string describing the data pair
    //! @return formatted string describing the data pair
    std::string DataPairwise::GetIdentification() const
    {
      // true if either of the locators is not defined
      if( !IsSet())
      {
        return std::string();
      }

      std::stringstream write;
      write << First()->GetIdentification();
      io::Serialize::Write( std::string( "<=>"), write, 1);
      write << Second()->GetIdentification();
      // return identification string
      return write.str();
    }

    //! @brief reads formatted string describing the locator
    //! @return formatted string describing the locator
    std::istream &DataPairwise::ReadIdentification( std::istream &ISTREAM)
    {
      m_First->ReadIdentification( ISTREAM);
      BCL_MessageDbg( "read in " + m_First->GetIdentification());
      std::string separator;
      io::Serialize::Read( separator, ISTREAM);
      BCL_MessageDbg( "read in " + separator);
      m_Second->ReadIdentification( ISTREAM);
      BCL_MessageDbg( "read in " + m_Second->GetIdentification());
      return ISTREAM;
    }

    //! @brief gives the first point in the data pair
    //! @return the first point in the data pair
    const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &DataPairwise::First() const
    {
      // true if this is not set
      if( !IsSet())
      {
        static const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> s_undefined;
        return s_undefined;
      }

      // return the first locator
      return m_First;
    }

    //! @brief gives the second point in the data pair
    //! @return the second point in the data pair
    const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &DataPairwise::Second() const
    {
      // true if this is not set
      if( !IsSet())
      {
        static util::ShPtr< assemble::LocatorAtomCoordinatesInterface> s_undefined;
        return s_undefined;
      }

      // return the second locator
      return m_Second;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief sets the two data points involved in this data pair
    //! @param ATOM_A the  first data point that will be set
    //! @param ATOM_B the second data point that will be set
    void DataPairwise::Set
    (
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &ATOM_A,
      const util::ShPtr< assemble::LocatorAtomCoordinatesInterface> &ATOM_B
    )
    {
      if( assemble::LocatorAtomCoordinatesInterface::PtrLessThan()( ATOM_A, ATOM_B))
      {
        m_First = ATOM_A;
        m_Second = ATOM_B;
      }
      else
      {
        m_First = ATOM_B;
        m_Second = ATOM_A;
      }
    }

    //! @brief indicates whether or not the two data points have been set
    //! @return boolean true if the data pair has been successfully set - false otherwise
    bool DataPairwise::IsSet() const
    {
      return m_First.IsDefined() && m_Second.IsDefined();
    }

    //! @brief calculates the euclidian distance indicated by a protein model
    //! @param MODEL the model from which the distance will be calculated
    //! @return double which is the distance between the two coordinates indicated by this data pair
    double DataPairwise::EuclidianDistance( const assemble::ProteinModel &MODEL) const
    {
      // locate the coordinates for the two locators
      const linal::Vector3D coords_a( m_First->Locate( MODEL));
      const linal::Vector3D coords_b( m_Second->Locate( MODEL));

      // true if either of the coordinates are not defined
      if( !coords_a.IsDefined() || !coords_b.IsDefined())
      {
        // return undefined double
        return util::GetUndefinedDouble();
      }

      // holds the distance between the two coordinate points
      const double distance( linal::Distance( coords_a, coords_b));

      // return the calculated distance
      return distance;
    }

    //! @brief calculates statistics for the distance indicated by this across an ensemble of models
    //! @param ENSEMBLE the ensemble over which distances and statistics will be calculated
    //! @return pair of RunningAverageSD< double> and RunningMinMax< double> indicated distance statistics over ENSEMBLE
    storage::Pair< math::RunningAverageSD< double>, math::RunningMinMax< double> >
    DataPairwise::EuclidianDistance( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // will keep track of the min max mean and stddev statistics
      storage::Pair< math::RunningAverageSD< double>, math::RunningMinMax< double> > stats;

      // iterate through the models of the ensemble
      for
      (
        assemble::ProteinEnsemble::const_iterator ensemble_itr( ENSEMBLE.Begin()), ensemble_itr_end( ENSEMBLE.End());
        ensemble_itr != ensemble_itr_end;
        ++ensemble_itr
      )
      {
        // locate the coordinates for the two locators
        linal::Vector3D coords_a( m_First->Locate( **ensemble_itr));
        linal::Vector3D coords_b( m_Second->Locate( **ensemble_itr));

        // holds the distance between the two coordinate points
        const double distance( linal::Distance( coords_a, coords_b));

        // true if either of the coordinates is not defined or the distance is not defined
        if( !coords_a.IsDefined() || !coords_b.IsDefined() || !util::IsDefined( distance))
        {
          // go to next model in ensemble
          continue;
        }

        // add the distance two the two statistic objects
        stats.First() += distance;
        stats.Second() += distance;
      }

      // return the statistics objects
      return stats;
    }

    //! @brief calculates the sequence separation between the two data points in this data pair
    //! @return size_t which is the sequence separation between the two data points in this data pair
    size_t DataPairwise::SequenceSeparation() const
    {
      return CalculateSequenceSeparation
      (
        First()->GetChainID(), First()->GetSeqID(), Second()->GetChainID(), Second()->GetSeqID()
      );
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataPairwise::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_First, ISTREAM);
      io::Serialize::Read( m_Second, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DataPairwise::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_First, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Second, OSTREAM, INDENT) << '\n';

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculates sequence separation between two residues indicated by chain and seq ids
    //! @param CHAIN_ID_A the chain id of the first residue
    //! @param SEQ_ID_A the seq id of the first residue
    //! @param CHAIN_ID_B the chain id of the second residue
    //! @param SEQ_ID_B the seq id of the second residue
    size_t DataPairwise::CalculateSequenceSeparation
    (
      const char CHAIN_ID_A, const int SEQ_ID_A, const char CHAIN_ID_B, const int SEQ_ID_B
    )
    {
      if( CHAIN_ID_A != CHAIN_ID_B)
      {
        return util::GetUndefinedSize_t();
      }

      return math::Absolute( SEQ_ID_A - SEQ_ID_B);
    }

    //! @brief less than operator for comparing two DataPairwise
    //! @param LHS the first DataPairwise which will be compared against the second DataPairwise
    //! @param RHS the second DataPairwise which will be compared against the first DataPairwise
    //! @return boolean true if LHS is less than RHS - false otherwise
    bool operator <( const DataPairwise &LHS, const DataPairwise &RHS)
    {
      // true if LHS is not set - cannot be less than RHS
      if( !LHS.IsSet())
      {
        BCL_MessageDbg( "LHS " + LHS.GetIdentification() + " RHS " + RHS.GetIdentification());
        return false;
      }

      // true if RHS is not set but LHS is set - LHS is less than RHS
      if( LHS.IsSet() && !RHS.IsSet())
      {
        BCL_MessageDbg( "LHS " + LHS.GetIdentification() + " RHS " + RHS.GetIdentification());
        return true;
      }

      // make vector nd 2 of locator atom from the LHS argument
      const storage::VectorND< 2, assemble::LocatorAtom> lhs
      (
        assemble::LocatorAtom( LHS.First()->GetChainID(), LHS.First()->GetSeqID(), LHS.First()->GetAtomType()),
        assemble::LocatorAtom( LHS.Second()->GetChainID(), LHS.Second()->GetSeqID(), LHS.Second()->GetAtomType())
      );

      // make vector nd 2 of locator atom from the RHS argument
      const storage::VectorND< 2, assemble::LocatorAtom> rhs
      (
        assemble::LocatorAtom( RHS.First()->GetChainID(), RHS.First()->GetSeqID(), RHS.First()->GetAtomType()),
        assemble::LocatorAtom( RHS.Second()->GetChainID(), RHS.Second()->GetSeqID(), RHS.Second()->GetAtomType())
      );

      BCL_MessageDbg( "LHS " + LHS.GetIdentification() + " RHS " + RHS.GetIdentification());
      // use the less than operator of the locator atom and the vector nd
      return lhs < rhs;
    }

    //! @brief writes pymol script formatted data to a stream that can show distances in pymol
    //! @param OSTREAM the stream which will write the data
    //! @param DISTANCES the data that will be written along with their distance and the desired color of their line
    //! @return ostream the stream that wrote the data distance information
    std::ostream &ShowDistancesInPymol
    (
      std::ostream &OSTREAM,
      const storage::List< storage::Triplet< DataPairwise, double, linal::Vector3D> > &DISTANCES
    )
    {
      // iterate through the distances
      for
      (
        storage::List< storage::Triplet< DataPairwise, double, linal::Vector3D> >::const_iterator
          distance_itr( DISTANCES.Begin()), distance_itr_end( DISTANCES.End());
        distance_itr != distance_itr_end;
        ++distance_itr
      )
      {
        // get critical data
        const std::string first_pymol_name( distance_itr->First().First()->GetPymolName());
        const std::string second_pymol_name( distance_itr->First().Second()->GetPymolName());
        const std::string dist_name( first_pymol_name + "_" + second_pymol_name);
        const std::string color_name( "color_" + dist_name);
        const char chain_a( distance_itr->First().First()->GetChainID());
        const int  resi_a(  distance_itr->First().First()->GetSeqID());
        const char chain_b( distance_itr->First().Second()->GetChainID());
        const int  resi_b(  distance_itr->First().Second()->GetSeqID());
        std::string atom_a( distance_itr->First().First()->GetAtomType());
        std::string atom_b( distance_itr->First().Second()->GetAtomType());

        // no hydrogens so replace with CA for visualization purposes
        const util::StringReplacement replacer
        (
          util::StringReplacement::e_Any, biol::GetAtomTypes().HA2, biol::GetAtomTypes().CA
        );
        replacer.ReplaceEachIn( atom_a);
        replacer.ReplaceEachIn( atom_b);

        // write the distance selection
        OSTREAM << "distance " << dist_name << ", (chain " << chain_a << " and resi " << resi_a << " and name " << atom_a
        << "),(chain " << chain_b << " and resi " << resi_b << " and name " << atom_b << ")" << '\n';

        // define the color for the current line
        const linal::Vector3D color( distance_itr->Third());
        OSTREAM << "set_color " << color_name << ", [" << color.X() << "," << color.Y() << "," << color.Z() << "]"
        << '\n';

        // set the color of the dashes
        OSTREAM << " set dash_color, " << color_name << ", " << dist_name << '\n';

        // set the gap between dashes
        static const double s_dash_gap( 0);
        OSTREAM << " set dash_gap, " << s_dash_gap << ", " << dist_name << '\n';

        // set the dash width
        static const double s_dash_width( 5);
        OSTREAM << " set dash_width, " << s_dash_width << ", " << dist_name << '\n';
      }

      return OSTREAM;
    }

  } // namespace restraint

} // namespace bcl
