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
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "descriptor/bcl_descriptor_molecule_entropy_qha.h"
#include "linal/bcl_linal_principal_component_analysis.h"
#include "math/bcl_math_running_average.h"
#include "quality/bcl_quality_rmsd.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoleculeEntropyQHA::MoleculeEntropyQHA()
    {}

    //! @brief
    //! @brief SAMPLER sample molecule conformations
    MoleculeEntropyQHA::MoleculeEntropyQHA
    (
      const chemistry::SampleConformations &SAMPLER
    ) :
      m_SampleConfs( SAMPLER)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeEntropyQHA
    MoleculeEntropyQHA *MoleculeEntropyQHA::Clone() const
    {
      return new MoleculeEntropyQHA( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeEntropyQHA::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeEntropyQHA::GetAlias() const
    {
      static const std::string s_name( "EntropyQHA");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeEntropyQHA::Calculate( linal::VectorReference< float> &STORAGE)
    {
      // initialize ensembles
      chemistry::FragmentEnsemble local_ensemble, global_ensemble;
      util::SiPtr< const chemistry::ConformationInterface> reference_mol( this->GetCurrentObject());
      chemistry::FragmentComplete mol( *reference_mol);
      util::SiPtrVector< const linal::Vector3D> coords( mol.GetAtomCoordinates());

      // get local ensemble
      if( m_LocalEnsembleFilename.size() <= 0)
      {
        m_SampleConfs.SetSamplingPreferences( false, false, true, false);
        local_ensemble = m_SampleConfs( mol).First();
      }

      // Get global ensemble
      if( m_GlobalEnsembleFilename.size() <= 0)
      {
        m_SampleConfs.SetSamplingPreferences( true, true, true, false);
        global_ensemble = m_SampleConfs( mol).First();
      }

      // Realign local ensemble
      for( auto itr( local_ensemble.Begin()), itr_end( local_ensemble.End()); itr != itr_end; ++itr)
      {
        auto rotamer_coords( itr->GetAtomCoordinates());
        auto transform( quality::RMSD::SuperimposeCoordinates( coords, rotamer_coords));
        itr->Transform( transform);
      }

      // Realign global ensemble
      for( auto itr( global_ensemble.Begin()), itr_end( global_ensemble.End()); itr != itr_end; ++itr)
      {
        auto rotamer_coords( itr->GetAtomCoordinates());
        auto transform( quality::RMSD::SuperimposeCoordinates( coords, rotamer_coords));
        itr->Transform( transform);
      }

      // remove hydrogen atoms
      for
      (
          auto local_itr( local_ensemble.Begin()), local_itr_end( local_ensemble.End());
          local_itr != local_itr_end;
          ++local_itr
      )
      {
        local_itr->RemoveH();
      }
      for
      (
          auto global_itr( global_ensemble.Begin()), global_itr_end( global_ensemble.End());
          global_itr != global_itr_end;
          ++global_itr
      )
      {
        global_itr->RemoveH();
      }

      // Make some matrices to fill in later, rows = number of conformers, cols = number of atoms in molecule
      linal::Matrix< double> local_coords( local_ensemble.GetSize(), size_t( this->GetCurrentObject()->GetSize() * 3));
      linal::Matrix< double> global_coords( global_ensemble.GetSize(), size_t( this->GetCurrentObject()->GetSize() * 3));

      BCL_MessageVrb( "Size of local conformational ensemble: " + util::Format()( local_ensemble.GetSize()));
      BCL_MessageVrb( "Size of global conformational ensemble: " + util::Format()( global_ensemble.GetSize()));

      // Fill in coordinates of local conformational ensemble matrix
      math::RunningAverage< linal::Vector< double>> row_ave;
      size_t conformer_index( 0);
      for
      (
          chemistry::FragmentEnsemble::iterator local_itr( local_ensemble.Begin()), local_itr_end( local_ensemble.End());
          local_itr != local_itr_end;
          ++local_itr, ++conformer_index
      )
      {
        size_t coords_index( 0);
        for( auto coords_itr( local_itr->GetAtomsIterator()); coords_itr.NotAtEnd(); ++coords_itr)
        {
          // mass reweighting
          double mass( math::Sqrt( coords_itr->GetElementType()->GetProperty( chemistry::ElementTypeData::e_Mass)));
          local_coords[conformer_index][coords_index++] = coords_itr->GetPosition().X() * mass;
          local_coords[conformer_index][coords_index++] = coords_itr->GetPosition().Y() * mass;
          local_coords[conformer_index][coords_index++] = coords_itr->GetPosition().Z() * mass;
        }
        row_ave += local_coords.GetRow( conformer_index);
      }

      // subtract average row value from each row element
      for( size_t row_index( 0); row_index < local_coords.GetNumberRows(); ++row_index)
      {
        auto row( local_coords.GetRow( row_index));
        row -= row_ave.GetAverage();
      }

      // Output the row norm for debugging purposes
      math::RunningAverage< double> l_ave;
      l_ave.Reset();
      linal::Vector< double> averages( local_coords.GetNumberRows());
      for( size_t i( 0); i < local_coords.GetNumberRows(); ++i)
      {
        auto row( local_coords.GetRow( i));
        for( size_t j( 0); j < row.GetSize(); ++j)
        {
          l_ave += row( j);
        }
        averages( i) = l_ave.GetAverage();
      }

      // Local ensemble correction:
      // Molecules with essentially only one allowable set of coordinates within
      // the prescribed dihedral bins will have eigenvalues of all 0, resulting in nan
      // So, if the row norm is 0 then add some noise
      if( averages.Norm() <= 1e-16)
      {
        for( size_t i( 0); i < local_coords.GetNumberRows(); ++i)
        {
          auto row( local_coords.GetRow( i));
          for( size_t j( 0); j < row.GetSize(); ++j)
          {
            double noise( random::GetGlobalRandom().RandomGaussian( 1e-11, 1e-12));
            row( j) += noise;
          }
        }
      }

      // Fill in coordinates of global conformational ensemble matrix
      row_ave.Reset();
      conformer_index = 0;
      for
      (
          chemistry::FragmentEnsemble::iterator global_itr( global_ensemble.Begin()), global_itr_end( global_ensemble.End());
          global_itr != global_itr_end;
          ++global_itr, ++conformer_index
      )
      {
        size_t coords_index( 0);
        for( auto coords_itr( global_itr->GetAtomsIterator()); coords_itr.NotAtEnd(); ++coords_itr)
        {
          double mass( math::Sqrt( coords_itr->GetElementType()->GetProperty( chemistry::ElementTypeData::e_Mass)));
          global_coords[conformer_index][coords_index++] = coords_itr->GetPosition().X() * mass;
          global_coords[conformer_index][coords_index++] = coords_itr->GetPosition().Y() * mass;
          global_coords[conformer_index][coords_index++] = coords_itr->GetPosition().Z() * mass;
        }
        row_ave += global_coords.GetRow( conformer_index);
      }

      // subtract average row value from each row element
      for( size_t row_index( 0); row_index < global_coords.GetNumberRows(); ++row_index)
      {
        auto row( global_coords.GetRow( row_index));
        row -= row_ave.GetAverage();
      }

      // Perform PCA on coordinate matrices
      linal::PrincipalComponentAnalysis< double> pca_object;
      storage::Pair< linal::Matrix< double>, linal::Vector< double> > local_pca = pca_object.GetSortedEigenVectorsValues(local_coords);
      storage::Pair< linal::Matrix< double>, linal::Vector< double> > global_pca = pca_object.GetSortedEigenVectorsValues(global_coords);

      // Eigenvalue sums
      double global_pca_sum = global_pca.Second().Sum();
      double local_pca_sum = local_pca.Second().Sum();
      double ratio_pca_sums = -std::log( local_pca_sum / global_pca_sum);

      // If the log ratio of the eigenvalue sums is greater than 4.0 then reset local_pca values to achieve 4.0
      bool rescale( ratio_pca_sums > 4.0);

      // Define extra constants
      double temp( 298.15); // K
      double h_kT = math::g_Plank / ( math::g_Boltzmann * temp) / 1.0e-12; // s

      // Convert eigenvalues to QHA frequencies
      size_t n_ev( local_pca.Second().GetSize());
      linal::Vector< double> local_freq( n_ev);
      size_t local_freq_index( 0);
      for
      (
          auto local_pca_itr( local_pca.Second().Begin()), local_pca_itr_end( local_pca.Second().End());
          local_pca_itr != local_pca_itr_end;
          ++local_pca_itr, ++local_freq_index
      )
      {
        if( rescale)
        {
          *local_pca_itr = *local_pca_itr * std::exp( ratio_pca_sums) / std::exp( 4.0);
        }
        if( *local_pca_itr > double( 0.0))
        {
          local_freq( local_freq_index) = math::Sqrt( std::max( ( math::g_Boltzmann * temp) / *local_pca_itr, 0.0));
        }
        else
        {
          local_freq( local_freq_index) = double( 0.0);
        }
      }

      linal::Vector< double> global_freq( global_pca.Second().GetSize());
      size_t global_freq_index( 0);
      for
      (
          auto global_pca_itr( global_pca.Second().Begin()), global_pca_itr_end( global_pca.Second().End());
          global_pca_itr != global_pca_itr_end;
          ++global_pca_itr, ++global_freq_index
      )
      {
        if( *global_pca_itr > double( 0.0))
        {
          global_freq( global_freq_index) = math::Sqrt( std::max( ( math::g_Boltzmann * temp) / *global_pca_itr, 0.0));
        }
        else
        {
          global_freq( global_freq_index) = double( 0.0);
        }
      }

      // Compute the local ensemble conformational entropy
      double local_s( 0.0);
      for( size_t local_freq_index( 0); local_freq_index < local_freq.GetSize(); ++local_freq_index)
      {
        if( local_freq( local_freq_index) > 0.0)
        {
          local_s += ( ( h_kT * local_freq( local_freq_index)) / ( std::exp( h_kT * local_freq( local_freq_index)) - 1)) - std::log( 1 - std::exp( -h_kT * local_freq( local_freq_index)));
        }
      }

      // Compute the global ensemble conformational entropy
      double global_s( 0.0);
      for( size_t global_freq_index( 0); global_freq_index < global_freq.GetSize(); ++global_freq_index)
      {
        if( global_freq( global_freq_index) > 0.0)
        {
          global_s += ( ( h_kT * global_freq( global_freq_index)) / ( std::exp( h_kT * global_freq( global_freq_index)) - 1)) - std::log( 1 - std::exp( -h_kT * global_freq( global_freq_index)));
        }
      }

      // compute final entropies
      global_s = global_s * math::g_Gas;
      local_s = local_s * math::g_Gas;

      // Store pertinent values
      STORAGE( 0) = global_s;
      STORAGE( 1) = local_s;
      STORAGE( 2) = global_s - local_s; // units of J/(K*mol) --> to get kcal/mol, multiply by (0.000239006 * 298.15)
      STORAGE( 3) = -std::log( local_s / global_s);

      STORAGE( 4) = global_pca.Second().Sum();
      STORAGE( 5) = local_pca.Second().Sum();
      STORAGE( 6) = global_pca.Second().Sum() - local_pca.Second().Sum();
      STORAGE( 7) = std::min( -std::log( local_pca.Second().Sum() / global_pca.Second().Sum()), 4.0);

      //Fix nonsense
      if( !STORAGE.CreateSubVectorReference( 4, 0).IsDefined() || STORAGE( 3) < 0.0)
      {
        STORAGE.CreateSubVectorReference( 4, 0) = 0.0;
      }
      if( !util::IsDefined( STORAGE( 7)))
      {
        STORAGE.CreateSubVectorReference( 3, 4) = 0.0;
        STORAGE( 7) = 4.0;
      }
    } // end Calculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeEntropyQHA::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeEntropyQHA::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Sampling options for relative conformational entropy estimates. "
          "Descriptor output indices correspond to the global_s, local_s, global_s - local_s, and -ln(local_s/global_s). "
          "The first 4 indices are actual entropy estimates in J*K^-1*mol^-1, while the latter 4 indices are the PCA eigenvalue sums.");
      parameters.AddInitializer
      (
        "sampler",
        "sample configurational space across dihedral bins",
        io::Serialization::GetAgent( &m_SampleConfs),
        util::ObjectDataLabel
        ( "("
          "conformation_comparer=bcl::chemistry::ConformationComparisonInterface,tolerance=0.0,"
          "generate_3D=0,relative_random_dihedral_change_weight=0.0,cluster=false,"
          "max_iterations=10000, max_conformations=1000"
          ")"
        )
      );
      parameters.AddInitializer
      (
        "local_ensemble",
        "input a pre-generated local ensemble",
        io::Serialization::GetAgent( &m_LocalEnsembleFilename),
        ""
      );
      parameters.AddInitializer
      (
        "global_ensemble",
        "input a pre-generated global ensemble",
        io::Serialization::GetAgent( &m_GlobalEnsembleFilename),
        ""
      );

      return parameters;
    }
  } // namespace descriptor
} // namespace bcl
