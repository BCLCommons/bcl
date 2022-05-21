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
#include "restraint/bcl_restraint_sas_debye.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_ca_cb.h"
#include "biol/bcl_biol_atom_group_types.h"
#include "biol/bcl_biol_sasa_data.h"
#include "fold/bcl_fold_add_parabolic_loops.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_sum_function.h"
#include "storage/bcl_storage_triplet.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SasDebye::s_Instance
    (
      util::Enumerated< SasDebyeInterface>::AddInstance( new SasDebye())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Constructor
    //! @param LOOPS bool value to represent loops that are not present in the protein model
    //! @param USE_REGULA_FALSI_APPROXIMATION whether to determine the norm factor with regula falsi (true) or
    //!        pythagorean approximation (false)
    //! @param EXCLUDED_VOLUME_PARAMETER C1 tuning parameter for data fitting
    //! @param HYDRATION_SHELL_PARAMETER C2 tuning parameter for data fitting
    //! @param SIDE_CHAIN_APPROXIMATION - true to approximate side chains, false otherwise
    //! @param REDUCED_EXPERIMENTAL_DATA - small data set of experimental data based on Shannon Sampling
    SasDebye::SasDebye
    (
      const bool LOOPS,
      const bool USE_REGULA_FALSI_APPROXIMATION,
      double EXCLUDED_VOLUME_PARAMETER,
      double HYDRATION_SHELL_PARAMETER,
      const bool SIDE_CHAIN_APPROXIMATION,
      const bool USE_SANS,
      double DEUTERIUM_EXCHANGE_PARAMETER,
      const util::ShPtr< storage::Vector< SasScatteringPoint> > REDUCED_EXPERIMENTAL_DATA
    )
    :
      m_ShouldApproximateLoops( LOOPS),
      m_DetermineAnalyticNormFactor( USE_REGULA_FALSI_APPROXIMATION),
      m_ExcludedVolumeParameter( EXCLUDED_VOLUME_PARAMETER),
      m_HydrationShellParameter( HYDRATION_SHELL_PARAMETER),
      m_ShouldApproximateSideChains( SIDE_CHAIN_APPROXIMATION),
      m_UseSans( USE_SANS),
      m_DeuteriumExchangeParameter( DEUTERIUM_EXCHANGE_PARAMETER)
    {
      SasDebyeInterface::SetReducedExperimentalData( REDUCED_EXPERIMENTAL_DATA);
    }

    //! @brief Clone function
    //! @return pointer to new SasDebye function
    SasDebye *SasDebye::Clone() const
    {
      return new SasDebye( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SasDebye::GetAlias() const
    {
      static const std::string s_Name( "SasDebye");
      return s_Name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SasDebye::GetClassIdentifier() const
    {
      // Get BCL Standardized Class name
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Get the positions and form factor functions for all the atoms in the protein model
    //! @param MODEL the protein model of interest
    //! @return a vector containing pairs of position and form factor function
    storage::Vector
    <
      storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
    > SasDebye::GetAtomsAndFormFactors( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      m_CompleteAllAtoms.Reset();
      m_AtomicCoordinates.Reset();
      m_AtomSasa.Reset();

      BCL_MessageDbg( " Inside Get Atoms and Form Factors Function: ");

      // initialize vector to hold atom coordinates and Form Factor F(q) function
      storage::Vector
      <
        storage::Triplet
        <
          linal::Vector3D,
          util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >,
          double
        >
      > all_atoms;

      size_t all_atoms_size( all_atoms.GetSize());
      BCL_MessageDbg( "initial all_atoms size: " + util::Format()( all_atoms_size));

      size_t data_member_size( m_AtomicCoordinates.GetSize());
      BCL_MessageDbg( "initial Data member size: " + util::Format()( data_member_size) + "\n");

      // Get Structure Factors for Loop coordinate approximation.  The form factors for all atoms in the residue will
      // be summed together
      if( m_ShouldApproximateLoops)
      {
        BCL_MessageDbg( " Inside the Loop approximation section: ");

        // storage for loop coordinates, if m_ShouldApproximateLoops was given
        const storage::Vector< storage::Pair< biol::AAType, linal::Vector3D> >
          loop_coordinates_aa_types( fold::AddParabolicLoops( m_DetermineAnalyticNormFactor).GetLoopCoordinates( PROTEIN_MODEL));
        for
        (
          storage::Vector< storage::Pair< biol::AAType, linal::Vector3D> >::const_iterator
            itr_loop_coords_aa_types( loop_coordinates_aa_types.Begin()),
            itr_loop_coords_aa_types_end( loop_coordinates_aa_types.End());
          itr_loop_coords_aa_types != itr_loop_coords_aa_types_end;
          ++itr_loop_coords_aa_types
        )
        {
          all_atoms.PushBack
          (
            storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
            (
              itr_loop_coords_aa_types->Second(),
              itr_loop_coords_aa_types->First()->GetStructureFactor(),
              0.0
            )
          );

          // Expand the residue to its atom types
          const storage::Set< biol::AtomType> &atoms_in_residue( itr_loop_coords_aa_types->First()->GetAllowedHeavyAtomTypes());
          for
           (
             storage::Set< biol::AtomType>::const_iterator
               itr_atom_type( atoms_in_residue.Begin()),
               itr_atom_type_end( atoms_in_residue.End());
             itr_atom_type != itr_atom_type_end;
             ++itr_atom_type
           )
          {

            biol::AAType aa_type( itr_loop_coords_aa_types->First());
            biol::AtomType atom_type( *itr_atom_type);

            m_CompleteAllAtoms.PushBack( biol::GetAtomGroupTypes().GetTypeString( aa_type, atom_type));

            m_AtomicCoordinates.PushBack( itr_loop_coords_aa_types->Second());
            m_AtomSasa.PushBack( 0.0);
          }

        } // end loop estimation
      }

      all_atoms_size = all_atoms.GetSize();
      BCL_MessageDbg( "after loop approximation all_atoms size: " + util::Format()( all_atoms_size));

      data_member_size = m_AtomicCoordinates.GetSize();
      BCL_MessageDbg( "after loop approximation Data member size: " + util::Format()( data_member_size) + "\n");

      size_t atom_number( 0);
      size_t sse_count( 0);
      size_t residue_count( 0);

      // iterate over the chains
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end; ++chain_itr
      )
      {
        // get all the structured sses
        const util::SiPtrVector< const assemble::SSE> structured_sses
        (
          ( *chain_itr)->GetSSEs( storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND))
        );

        // if structured_sses are empty
        if( structured_sses.IsEmpty() && m_ShouldApproximateLoops)
        {
          // warn user and return empty vector
          BCL_Message( util::Message::e_Standard, "No structured SSEs found in protein model");
          continue;
        }

        // iterate over all the SSEs including coil (depending on min_sse_size parameter)
        for
        (
          storage::Set< util::ShPtr< assemble::SSE>, assemble::SSELessThanNoOverlap>::const_iterator
            sse_itr( ( *chain_itr)->GetData().Begin()), sse_itr_end( ( *chain_itr)->GetData().End());
          sse_itr != sse_itr_end; ++sse_itr
        )
        {

           BCL_MessageDbg( "sse_count: " + util::Format()( sse_count));

          // iterate over the Amino Acids in a given SSE
          for
          (
            util::ShPtrVector< biol::AABase>::const_iterator aa_itr( ( *sse_itr)->GetData().Begin()),
              aa_itr_end( ( *sse_itr)->GetData().End());
            aa_itr != aa_itr_end; ++aa_itr
          )
          {
            ++residue_count;
            std::string amino_acid( ( *aa_itr)->GetType()->GetName());
            //BCL_MessageStd( "Residue: " + util::Format()( amino_acid));

            // initialize sum function beta carbon and for missing atom structure factors if specified
            util::ShPtr< math::SumFunction< SasDataParameters, double> > cb_atom_structure_factors
            (
              new math::SumFunction< SasDataParameters, double>()
            );

            // get the atom types for this residue
            // Be sure to get the Heavy Atom Types ( all atoms except hydrogen).  This set contains all of the residues
            // That a given amino acid contains.  Not the atoms actually present in the input PDB file
            const storage::Set< biol::AtomType> &complete_atom_types( ( *aa_itr)->GetType()->GetAllowedHeavyAtomTypes());

            //BCL_MessageStd( "AtomTypes: " + util::Format()( complete_atom_types));

            // get the first side chain atom
            const biol::Atom &first_side_chain_atom( ( *aa_itr)->GetFirstSidechainAtom());

            // if you do not have experimental Sasa Data then atom_sasa will always be zero.
            double atom_sasa( 0.0);
            double beta_sasa( 0.0);

            // iterate over the atoms in a given residue
            for
            (
              storage::Set< biol::AtomType>::const_iterator atom_itr( complete_atom_types.Begin()),
                atom_itr_end( complete_atom_types.End());
              atom_itr != atom_itr_end; ++atom_itr
            )
            {
              atom_number++;

              // get the atom type
              const biol::AtomType atom_type_name( ( *aa_itr)->GetAtom( *atom_itr).GetType());

              // get the atom coordinates
              const linal::Vector3D &atom_coords( ( *aa_itr)->GetAtom( *atom_itr).GetCoordinates());

              // get the SASA data
              util::ShPtr< biol::SasaData> sp_sasa_data
              (
                PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Sasa)
              );

              if( atom_coords.IsDefined())
              {
                std::string atom( ( *atom_itr)->GetName());
              }

              // Using Solvent Excluded Surface Area
              // Make sure the atom coordinates are defined
              if( sp_sasa_data.IsDefined() && atom_coords.IsDefined())
              {
                // Get the pdbid for the atom type
                const size_t atom_id( ( *aa_itr)->GetAtom( *atom_itr).GetPdbID());
                BCL_MessageDbg( "atom_id: " + util::Format()( atom_id));

                atom_sasa = sp_sasa_data->GetData()[ atom_id - 1]->GetSolventAccessibleSurface();
                BCL_MessageDbg( "atom_sasa: " + util::Format()( atom_sasa));
              }

              // if the coordinates are defined and this is not CB
              if( atom_coords.IsDefined() && *atom_itr != first_side_chain_atom.GetType())
              {

                const std::string atomtypecheck( biol::GetAtomGroupTypes().GetTypeString( ( *aa_itr)->GetType(), *atom_itr));
                //BCL_MessageStd( "After map AtomTypes: " + util::Format()( atomtypecheck));

                // pushback into vector
                all_atoms.PushBack
                (
                  storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
                  (
                    atom_coords,
                    util::CloneToShPtr( *biol::GetAtomGroupTypes().GetType( ( *aa_itr)->GetType(), *atom_itr)),
                    atom_sasa
                  )
                );

                m_CompleteAllAtoms.PushBack
                (
                  biol::GetAtomGroupTypes().GetTypeString( ( *aa_itr)->GetType(), *atom_itr)
                );

                m_AtomicCoordinates.PushBack( atom_coords);

                ( sp_sasa_data.IsDefined()) ? m_AtomSasa.PushBack( atom_sasa / 100) : m_AtomSasa.PushBack( 0.0);

                all_atoms_size = all_atoms.GetSize();
                BCL_MessageDbg( "Line 323 all_atoms size: " + util::Format()( all_atoms_size));

                data_member_size = m_AtomicCoordinates.GetSize();
                BCL_MessageDbg( "Line 326 Data member size: " + util::Format()( data_member_size) + "\n");
              }
              else if( first_side_chain_atom.GetCoordinates().IsDefined())
              {
                // add cb atom to structure factor calculated
                if( *atom_itr == first_side_chain_atom.GetType())
                {
                  *cb_atom_structure_factors = *biol::GetAtomGroupTypes().GetType( ( *aa_itr)->GetType(), *atom_itr);

                  m_CompleteAllAtoms.PushBack
                  (
                    biol::GetAtomGroupTypes().GetTypeString( ( *aa_itr)->GetType(), *atom_itr)
                  );

                  m_AtomicCoordinates.PushBack( first_side_chain_atom.GetCoordinates());
                  linal::Vector3D tempcoords( first_side_chain_atom.GetCoordinates());

                  ( sp_sasa_data.IsDefined()) ? m_AtomSasa.PushBack( atom_sasa / 100) : m_AtomSasa.PushBack( 0.0);

                  // Store Atom Sasa at this iteration to insert into CB later
                  beta_sasa = atom_sasa;

                  all_atoms_size = all_atoms.GetSize();
                  BCL_MessageDbg( "Line 349 all_atoms size: " + util::Format()( all_atoms_size));

                  data_member_size = m_AtomicCoordinates.GetSize();
                  BCL_MessageDbg( "Line 352 Data member size: " + util::Format()( data_member_size) + "\n:");
                }

                // add missing residues to structure factor calculation for CB
                if( m_ShouldApproximateSideChains && *atom_itr != first_side_chain_atom.GetType())
                {
                  *cb_atom_structure_factors += *biol::GetAtomGroupTypes().GetType( ( *aa_itr)->GetType(), *atom_itr);

                  // This is the case where the atom coordinates are not defined and the form factors are summed on
                  // the CB position. In this situation the SASA value will be 0

                  m_CompleteAllAtoms.PushBack
                  (

                    biol::GetAtomGroupTypes().GetTypeString( ( *aa_itr)->GetType(), *atom_itr)
                  );

                  m_AtomicCoordinates.PushBack( first_side_chain_atom.GetCoordinates());
                  linal::Vector3D tempcoords( first_side_chain_atom.GetCoordinates());

                  m_AtomSasa.PushBack( 0.0);

                  all_atoms_size = all_atoms.GetSize();
                  BCL_MessageDbg( "Line 373 all_atoms size: " + util::Format()( all_atoms_size));

                  data_member_size = m_AtomicCoordinates.GetSize();
                  BCL_MessageDbg( "Line 376 Data member size: " + util::Format()( data_member_size));
                }
              }
            }

            // if cb is defined push back the cb point that will either be the CB atom alone, or the CB atom with
            // the other side chain atoms form factor contributions summed together
            if( !cb_atom_structure_factors->GetFunction().IsEmpty())
            {
              all_atoms.PushBack
              (
                storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
                (
                  first_side_chain_atom.GetCoordinates(),
                  cb_atom_structure_factors,
                  beta_sasa
                )
              );

              all_atoms_size = all_atoms.GetSize();
              BCL_MessageDbg
              (
                "Added Cb Atom to all_atoms: " + util::Format()( all_atoms_size) +
                "\n ------------------------------------------------------------------------- \n"
              );

            } // Atom Type Iteration
          } // AA iteration
        } // SSE iteration
      } // chain iteration

        all_atoms_size = all_atoms.GetSize();
        BCL_MessageDbg( "Line 400 all_atoms size: " + util::Format()( all_atoms_size));

        data_member_size = m_AtomicCoordinates.GetSize();
        BCL_MessageDbg
        ( "Line 403 Data member size: " + util::Format()( data_member_size) + "\n");

      return all_atoms;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a proteinModel as input and calculates an intensity using the debye formula
    //! @param PROTEIN_MODEL
    //! @return the intensity for a given q value
    SasExperimentalAndCalculatedData SasDebye::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {

//      //********************************** Block to circumvent SAXS Calculation *********************************
//
//      // Read Processed Data File and return result
//      util::ShPtr< restraint::SasExperimentalAndCalculatedData> sp_precalculated_data( new restraint::SasExperimentalAndCalculatedData);
//
//      io::IFStream read;
//      io::File::MustOpenIFStream( read, "scaled_y_axis.data");
//      sp_precalculated_data->ReadFromDataFile( read);
//      io::File::CloseClearFStream( read);
//
//      return *sp_precalculated_data;
//
//      //*************************************** End of Temp Block

      // get the experimental SAXS data
      util::ShPtr< SasScatteringData> sp_experimental_data;

      if( !this->GetReducedExperimentalData().IsDefined())
      {
        // get the experimental SAXS data
        sp_experimental_data = this->GetExperimentalData();

        // Verify you have experimental Saxs Data
        if( !sp_experimental_data.IsDefined())
        {
          // warn user and return empty data
          BCL_Message( util::Message::e_Critical, "No experimental SAXS data found, returning empty data");
          return SasExperimentalAndCalculatedData();
        }
      }
      else
      {
        SasScatteringData experimental_data;

        for
        (
          storage::Vector< SasScatteringPoint>::const_iterator
          exp_data_itr( this->GetReducedExperimentalData()->Begin()),
          exp_data_itr_end( this->GetReducedExperimentalData()->End());
          exp_data_itr != exp_data_itr_end;
          ++exp_data_itr
        )
        {
          experimental_data.PushBackScattering( *exp_data_itr);
        }

        // use the Clone to ShPtr to create a shared pointer to experimental data
        util::ShPtr< SasScatteringData> sp_reduced_data( util::CloneToShPtr( experimental_data));

        // set sp_experimental_data to the reduced data set
        sp_experimental_data = sp_reduced_data;
      }

      storage::Vector
      <
        storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
      > all_atoms( GetAtomsAndFormFactors( PROTEIN_MODEL));

      // create object to hold calculated data
      SasScatteringData calculated_data;
      calculated_data.AllocateScatteringMemory( sp_experimental_data->GetScatteringData().GetSize());

      const size_t number_of_atoms( all_atoms.GetSize());
      linal::Matrix< double> distance_matrix( number_of_atoms, number_of_atoms);
      size_t row( 0);

      // Compute the Distance Matrix

      for
      (
        storage::Vector
        <
          storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
        >::const_iterator atom_itr_a( all_atoms.Begin()), atom_itr_end( all_atoms.End());
        atom_itr_a != atom_itr_end; ++atom_itr_a, ++row
      )
      {
        size_t col( row + 1);
        for
        (
          storage::Vector
          <
            storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
          >::const_iterator atom_itr_b( atom_itr_a + 1);
          atom_itr_b != atom_itr_end; ++atom_itr_b, ++col
        )
        {
          distance_matrix( row, col) = linal::Distance( atom_itr_a->First(), atom_itr_b->First());
        }
      }

      // iterate over experimental data to get q-values
      for
      (
        storage::Vector< SasScatteringPoint>::const_iterator
          data_itr( sp_experimental_data->GetScatteringData().Begin()),
          data_itr_end( sp_experimental_data->GetScatteringData().End());
        data_itr != data_itr_end;
        ++data_itr
      )
      {
        // variable to hold q-value
        const double &q( data_itr->GetQvalue());

        // setup variable to hold intensity
        double intensity( 0.0);
        storage::Vector< double> residues;
        residues.AllocateMemory( number_of_atoms);

        for
        (
          storage::Vector
          <
            storage::Triplet< linal::Vector3D, util::ShPtr< math::FunctionInterfaceSerializable< SasDataParameters, double> >, double>
          >::const_iterator atom_itr_a( all_atoms.Begin()), atom_itr_end( all_atoms.End());
          atom_itr_a != atom_itr_end; ++atom_itr_a
        )
        {
          // restraint::SasDataParameters data( q, atom_itr_a->Third(), m_ExcludedVolumeParameter, m_HydrationShellParameter);
          SasDataParameters data
          (
            m_UseSans,
            q,
            atom_itr_a->Third(),
            m_ExcludedVolumeParameter,
            m_HydrationShellParameter,
            m_DeuteriumExchangeParameter
          );

          // calculate form factor for given q value for residue type
          residues.PushBack( atom_itr_a->Second()->operator()( data));
        }

        // iterate over pairs in the vector
        for( size_t row( 0); row < number_of_atoms; ++row)
        {
          // Get residue a for the debye calculation
          const double residue_a( residues( row));

          // calculate I using cb-cb distance (rij) and form factors from the amino acid types
          // I(q) = (sum i=1 to M)(sum j=1 to M) Fi(q)*Fj(q) * sin(q*rij)/(q*rij)
          // which can be rewritten as
          // I(q) = (sum i=1 to M) Fi(q)^2 +  2 * Fi(q) * (sum j=i+1 to M) Fj(q) * sin(q*rij)/(q*rij)
          // here we will call (sum j=i+1 to M) Fj(q) * sin(q*rij)/(q*rij) inner sum
          double inner_sum( 0.0);

          // iterate over upper triangle of atom-atom distance pairs
          for( size_t col( row + 1); col < number_of_atoms; ++col)
          {

            // const double cb_distance( linal::Distance( atom_itr_a->First(), atom_itr_b->First()));

            // calculate form factor for given q value for residue type
            //const double residue_b( atom_itr_b->Second()->operator()( q));
            const double residue_b( residues( col));
            const double &current_distance = distance_matrix( row, col);

            // calculate x = q * current_distance
            const double x( q * current_distance);

            // Sum the intensity over the protein
            // There is a nuance at point q = 0.  The limit of the function sin(x)/x is 1 as x approaches 0.
            // Therefore, at q==0 we only sum by residue_b.  The second term goes to 1.
            if( q == 0.0)
            {
              inner_sum += residue_b;
            }
            else
            {
              inner_sum += residue_b * ( sin( x)) / x;

            }
          }
          intensity += 2.0 * residue_a * inner_sum + math::Sqr( residue_a);
        }

        // push back into storage::Vector< storage::VectorND< 3, double> >
        // first -> q, second -> I, third -> error
        calculated_data.PushBackScattering( SasScatteringPoint( q, intensity, 0.0));

      }
      SasExperimentalAndCalculatedData saxs_data( *sp_experimental_data, calculated_data);

      return saxs_data;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SasDebye::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "performs saxs debye calculation"
      );

      parameters.AddInitializer
      (
        "consider loops",
        "should loops be considered",
        io::Serialization::GetAgent( &m_ShouldApproximateLoops),
        "0"
      );
      parameters.AddInitializer
      (
        "analytic",
        "whether to determine the norm factor with regula falsi (1) or pythagorean approximation (0)",
        io::Serialization::GetAgent( &m_DetermineAnalyticNormFactor),
        "0"
      );

      parameters.AddInitializer
      (
        "excluded volume",
        "tuning parameter for excluded solvent volume",
        io::Serialization::GetAgent( &m_ExcludedVolumeParameter),
        "1.0"
      );

      parameters.AddInitializer
      (
        "hydration shell",
        "tuning parameter for hydration shell",
        io::Serialization::GetAgent( &m_HydrationShellParameter),
        "0.0"
      );

      parameters.AddOptionalInitializer
      (
        "approximate_sidechains",
        "sum up form factor contribution on cb position",
        io::Serialization::GetAgent( &m_ShouldApproximateSideChains)
      );

      parameters.AddOptionalInitializer
      (
        "use_sans",
        "use sans implementation of debye formula",
        io::Serialization::GetAgent( &m_UseSans)
      );

      parameters.AddInitializer
      (
        "deuterium_percentage",
        "percentage of deuterium comprising solvent",
        io::Serialization::GetAgent( &m_DeuteriumExchangeParameter),
        "0.0"
      );

      return parameters;
    }

  } // namespace restraint
} // namespace bcl
