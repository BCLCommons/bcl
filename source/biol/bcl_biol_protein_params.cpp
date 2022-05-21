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
#include "biol/bcl_biol_protein_params.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_protein_charge.h"
#include "opti/bcl_opti_approximator_root_regula_falsi.h"
#include "opti/bcl_opti_criterion_convergence_argument.h"
#include "storage/bcl_storage_table.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinParams::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinParams())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinParams::ProteinParams()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProteinParams
    ProteinParams *ProteinParams::Clone() const
    {
      return new ProteinParams( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinParams::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate protein parameters
    //! @param SEQUENCE amino acid sequence
    //! @return table containing all information
    storage::Table< double> ProteinParams::operator()( const AASequence &SEQUENCE) const
    {
      static const double s_weight_water
      (
        2 * chemistry::GetElementTypes().e_Hydrogen->GetProperty( chemistry::ElementTypeData::e_Mass) +
        1 * chemistry::GetElementTypes().e_Oxygen->GetProperty( chemistry::ElementTypeData::e_Mass)
      );

      // result
      storage::Table< double> result( storage::TableHeader( storage::Vector< std::string>::Create( "value", "percentage")));

      // count the number of amino acids of each type
      const storage::Map< AAType, size_t> number_aas( CountAAs( SEQUENCE));

      // add the aa count
      AddAACountToResultTable( number_aas, result);

      // molecular weight
      storage::Row< double> &mol_weight_row( result.InsertRow( "molecular_weight[Da]"));
      mol_weight_row( 0) = CalcualteMolecularWeight( number_aas) + s_weight_water;
      mol_weight_row( 1) = 100;

      // extinction coefficient
      CalcualteExtinctionCoefficient( number_aas, result);

      if( SEQUENCE.GetSize() < 1)
      {
        return result;
      }

      // pI
      for
      (
        AATypeData::PropertyTypeEnum pk_property( ProteinCharge::s_FirstpKAAProperty);
        pk_property <= ProteinCharge::s_LastpKAAProperty;
        ++pk_property
      )
      {
        storage::Row< double> &pi_row( result.InsertRow( "pI_" + pk_property.GetString()));
        pi_row( 0) = CalculatePI( number_aas, SEQUENCE.GetFirstAA()->GetType(), SEQUENCE.GetLastAA()->GetType(), pk_property);
        pi_row( 1) = 100;
      }

      // end
      return result;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinParams::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinParams::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate the molecular weight from the number of amino acids
    //! @param AA_COUNT a count for each amino acid type
    //! @return molecular weight in Dalton [Da]
    double ProteinParams::CalcualteMolecularWeight( const storage::Map< AAType, size_t> &AA_COUNT)
    {
      double mol_weight( 0.0);

      for
      (
        storage::Map< AAType, size_t>::const_iterator itr( AA_COUNT.Begin()), itr_end( AA_COUNT.End());
        itr != itr_end;
        ++itr
      )
      {
        mol_weight += itr->second * itr->first->GetAAProperty( AATypeData::e_Mass);
      }

      return mol_weight;
    }

    //! @brief calculate the pI value (the pH where the protein is not charged)
    //! @param AA_COUNT a count for each amino acid type
    //! @param N_TERM_TYPE type of n terminal amino acid
    //! @param C_TERM_TYPE type of c terminal amino acid
    //! @param PK_PROPERTY the pk value scale to use
    //! @return the PI value - the pH at which the protein would have no net-charge
    double ProteinParams::CalculatePI
    (
      const storage::Map< AAType, size_t> &AA_COUNT,
      const AAType &N_TERM_AATYPE,
      const AAType &C_TERM_AATYPE,
      const AATypeData::PropertyType &PK_PROPERTY
    )
    {
      // borders of the interval - ph 0 to 14
      const double border_left( 0);
      const double border_right( 14);

      // tolerance for the approximation result
      const double root_tolerance( 0.00001);

      // create termination criteria
      opti::CriterionConvergenceArgument< double, double> criterion( 1, root_tolerance);

      util::ShPtr< ProteinCharge> protein_charge_function( new ProteinCharge( AA_COUNT, N_TERM_AATYPE, C_TERM_AATYPE));
      protein_charge_function->SetPKProperty( PK_PROPERTY);

      // create the approximator
      opti::ApproximatorRootRegulaFalsi< double, double> approximator
      (
        *protein_charge_function, criterion, border_left, border_right
      );

      // approximation
      approximator.Approximate();
      const util::ShPtr< storage::Pair< double, double> > sp_result( approximator.GetTracker().GetBest());

      return sp_result->First();
    }

    //! @brief calculate the extinction coefficient from the number of amino acids in units of  M-1 cm-1, at 280 nm measured in water
    //! @param AA_COUNT a count for each amino acid type
    //! @param RESULT result table
    void ProteinParams::CalcualteExtinctionCoefficient
    (
      const storage::Map< AAType, size_t> &AA_COUNT, storage::Table< double> &RESULT
    )
    {
      const storage::Map< AAType, double> coefficient_map( ExtinctionCoefficientMap());

      double coefficient( 0);
      double coefficient_cystines( 0);
      for
      (
        storage::Map< AAType, double>::const_iterator itr( coefficient_map.Begin()), itr_end( coefficient_map.End());
        itr != itr_end;
        ++itr
      )
      {
        if( AA_COUNT.Has( itr->first))
        {
          if( itr->first == GetAATypes().CYS)
          {
            // only even number of cysteins considered, since only peptide bond does contribute
            coefficient_cystines += ( AA_COUNT.Find( itr->first)->second / 2) * ( 2 * itr->second);
            continue;
          }
          else
          {
            coefficient += AA_COUNT.Find( itr->first)->second * itr->second;
          }
        }
      }

      RESULT.InsertRow( "ExtinctionCoefficient[M-1*cm-1]")( 0) = coefficient;
      RESULT.InsertRow( "ExtinctionCoefficientCystines[M-1*cm-1]")( 0) = coefficient + coefficient_cystines;
    }

    //! @brief count the number of each amino acid in the given sequence
    //! @param SEQUENCE amino acid sequence
    //! @return the number of each amino acid type
    storage::Map< AAType, size_t> ProteinParams::CountAAs( const AASequence &SEQUENCE)
    {
      // count the number of amino acids of each type
      storage::Map< AAType, size_t> number_aas;

      for( AASequence::const_iterator itr( SEQUENCE.Begin()), itr_end( SEQUENCE.End()); itr != itr_end; ++itr)
      {
        ++number_aas[ ( *itr)->GetType()];
      }

      // end
      return number_aas;
    }

    //! @brief extinction coefficient map at 280 nm in water
    //! @return map of amino acid extinction coefficients
    storage::Map< AAType, double> ProteinParams::ExtinctionCoefficientMap()
    {
      storage::Map< AAType, double> map;

      map[ GetAATypes().TRP] = 5500.0;
      map[ GetAATypes().TYR] = 1490.0;
      map[ GetAATypes().CYS] =   62.5; // 125 for each cystine - which are two sulifde-bonded cysteins

      return map;
    }

    //! @brief add the aa count to the result table
    //! @param AA_COUNT map of aatype and their counts
    //! @param RESULT result table
    void ProteinParams::AddAACountToResultTable
    (
      const storage::Map< AAType, size_t> &AA_COUNT,
      storage::Table< double> &RESULT
    )
    {
      storage::Map< AAType, size_t> aa_count( AA_COUNT);
      size_t number_aas( 0);

      // iterate through all natural aa types
      for
      (
        AATypes::const_iterator itr( GetAATypes().Begin()), itr_end( GetAATypes().VAL.GetIterator() + 1);
        itr != itr_end;
        ++itr
      )
      {
        const AAType &current_type( *itr);
        const std::string name( current_type->GetThreeLetterCode() + '_' + current_type->GetOneLetterCode());
        storage::Row< double> &new_row( RESULT.InsertRow( name));
        new_row( 0) = aa_count[ current_type];
        number_aas += new_row( 0);
        aa_count.Erase( current_type);
      }

      // add the count of non-natural aa types
      for
      (
        storage::Map< AAType, size_t>::const_iterator itr( aa_count.Begin()), itr_end( aa_count.End());
        itr != itr_end;
        ++itr
      )
      {
        const AAType &current_type( itr->first);
        const std::string name( current_type->GetThreeLetterCode() + '_' + current_type->GetOneLetterCode());
        storage::Row< double> &new_row( RESULT.InsertRow( name));
        new_row( 0) = itr->second;
        number_aas += itr->second;
      }

      // total
      storage::Row< double> &new_row( RESULT.InsertRow( "number_aas"));
      new_row( 0) = number_aas;

      // percentage
      for
      (
        storage::Table< double>::iterator itr( RESULT.Begin()), itr_end( RESULT.End());
        itr != itr_end && itr->First() != "number_aas";
        ++itr
      )
      {
        itr->Second()( 1) = size_t( ( itr->Second()( 0) * 1000) / number_aas) / 10.0;
      }

      // neg charged residues
      {
        size_t neg_count( 0);
        if( AA_COUNT.Has( GetAATypes().ASP))
        {
          neg_count += AA_COUNT.Find( GetAATypes().ASP)->second;
        }
        if( AA_COUNT.Has( GetAATypes().GLU))
        {
          neg_count += AA_COUNT.Find( GetAATypes().GLU)->second;
        }

        storage::Row< double> &new_row( RESULT.InsertRow( "number_negative_aas"));
        new_row( 0) = neg_count;
        new_row( 1) = size_t( ( neg_count * 1000) / number_aas) / 10.0;
      }

      // pos charged residues
      {
        size_t pos_count( 0);
        if( AA_COUNT.Has( GetAATypes().ARG))
        {
          pos_count += AA_COUNT.Find( GetAATypes().ARG)->second;
        }
        if( AA_COUNT.Has( GetAATypes().LYS))
        {
          pos_count += AA_COUNT.Find( GetAATypes().LYS)->second;
        }

        storage::Row< double> &new_row( RESULT.InsertRow( "number_positive_aas"));
        new_row( 0) = pos_count;
        new_row( 1) = size_t( ( pos_count * 1000) / number_aas) / 10.0;
      }
    }

  } // namespace biol
} // namespace bcl
