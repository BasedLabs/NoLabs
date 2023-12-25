<template>
  <div class="text-center" v-if="this.$route.path === '/'">
    <div class="container-fluid">
      <div class="row">
        <!-- Amino Acid Lab -->
        <div class="lab-container col-md-4 col-xs-1">
          <div class="protein-container" ref="aminoAcidProtein"></div>
          <ul class="lab-features">
            <li>Folding</li>
            <li>Gene Ontology</li>
            <li>Solubility</li>
            <li>Localization</li>
          </ul>
          <RouterLink type="button" class="btn btn-primary btn-md" to="/amino-acid-lab">Amino acid lab</RouterLink>
        </div>

        <!-- Drug Discovery Lab -->
        <div class="lab-container col-md-4 col-xs-1">
          <div class="protein-container" ref="drugDiscoveryProtein"></div>
          <ul class="lab-features">
            <li>Protein-Ligand Binding</li>
            <li class="upcoming-feature">Metabolic processes affected</li>
          </ul>
          <RouterLink type="button" class="btn btn-primary btn-md" to="/drug-target-lab">Drug discovery lab</RouterLink>
        </div>

        <!-- Conformations lab Lab -->
        <div class="lab-container col-md-4 col-xs-1">
          <div class="protein-container" ref="conformationsProtein"></div>
          <ul class="lab-features">
            <li class="upcoming-feature">Protein conformations simulation</li>
          </ul>
          <RouterLink type="button" class="btn btn-primary btn-md" to="/conformations">Conformations lab (pre-alpha)
          </RouterLink>
        </div>
      </div>
      <div class="row">
        <!-- Protein viewer-->
        <div class="lab-container col-md-4 col-xs-1">
          <div class="protein-container"><img class="protein-container-gif" src="exampleProteins/shadi_3P.gif"></div>
          <RouterLink type="button" class="btn btn-primary btn-md" to="/protein-viewer">Simple .pdb viewer</RouterLink>
        </div>

        <!-- Protein motif Lab -->
        <div class="lab-container col-md-4 col-xs-1">
          <div class="protein-container" ref="proteinDesign"></div>
          <RouterLink type="button" class="btn btn-primary btn-md" to="/protein-design">Protein design lab</RouterLink>
        </div>
      </div>
    </div>
  </div>
  <RouterView />
</template>


<script>
import * as $3Dmol from '3dmol';

export default {
  watch: {
    '$route'(to, from) {
      // This will be called also when the route changes and the component is reused
      if (to.path === '/') {
        this.initProteins();
      }
    }
  },
  data() {
    return {
      dropdowns: {
        aminoAcid: false,
        drugDiscovery: false
      }
    };
  },
  methods: {
    initProteins() {
      this.loadProteinModel(this.$refs.aminoAcidProtein, 'exampleProteins/protein.pdb', true, true);
      this.loadProteinModel(this.$refs.drugDiscoveryProtein, 'exampleProteins/protein_ligand.pdb', false, true);
      this.loadProteinModel(this.$refs.conformationsProtein, 'exampleProteins/conformations.pdb', true, true);
      this.loadProteinModel(this.$refs.proteinDesign, 'exampleProteins/proteinDesign.pdb', true, false);
    },
    loadProteinModel(element, pdbPath, isCartoon, rotate) {
      let viewer = $3Dmol.createViewer(element, {
        backgroundColor: 'black'
      });

      fetch(pdbPath)
        .then(response => response.text())
        .then(pdbData => {
          viewer.addModel(pdbData, 'pdb');

          if (isCartoon) {
            viewer.setStyle({}, { cartoon: { color: 'spectrum' } });
          }
          else {
            viewer.setStyle({}, { cartoon: { color: 'white' } });
            // Style for the ligands
            // You can customize this style as needed
            viewer.setStyle({ hetflag: true }, { stick: { radius: 0.6, colorscheme: 'greenCarbon' } });
          }
          viewer.zoomTo();
          if(rotate){
            this.rotateModel(viewer);
          }
          viewer.render();
        });
    },
    rotateModel(viewer) {
      function rotate() {
        viewer.rotate(0.3, { x: 0, y: 1, z: 0 }); // Adjust rotation axis and speed as needed
        viewer.render();
        requestAnimationFrame(rotate);
      }
      rotate();
    },
    toggleDropdown(dropdown) {
      this.dropdowns[dropdown] = !this.dropdowns[dropdown];
    }
  }
};
</script>

<style>
.btn {
  /* Your existing styles for buttons */
  border-radius: 10px;
  /* Adjust the radius as needed */
}

.lab-features li {
  margin-bottom: 5px;
  /* Space between list items */
  position: relative;
  /* Needed for positioning the pseudo-element */
  padding-left: 25px;
  /* Space for the checkmark */
  text-align: left;
}

.lab-features li::before {
  content: '✅';
  /* The checkmark */
  position: absolute;
  left: 0;
  /* Aligns the checkmark to the left */
  top: 0;
  /* Aligns the checkmark to the top */
}

.lab-features .upcoming-feature::before {
  content: '🔜';
  /* Change the marker to "🔜" */
}

.toggle-text {
  cursor: pointer;
  color: #007bff;
  margin-top: 5px;
  margin-bottom: 5px;
  text-decoration: underline;
}

.dropdown-content {
  /* Your styles for the dropdown content */
  margin-top: 10px;
  border: 1px solid #ddd;
  padding: 10px;
  border-radius: 5px;
  box-shadow: 0 2px 5px rgba(0, 0, 0, 0.2);
}
</style>