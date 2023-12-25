<template>
  <div class="text-center" v-if="this.$route.path === '/'">
    <div class="container-fluid row justify-content-center">
      <!-- Amino Acid Lab -->
      <div class="col-md-3 d-flex flex-column align-items-center">
        <div class="mx-auto protein-container" ref="aminoAcidProtein"></div>
        <div class="lab-features">
          <span class="tag">Folding</span>
          <span class="tag">Gene Ontology</span>
          <span class="tag">Solubility</span>
          <span class="tag">Localization</span>
          <span class="upcoming-feature">Protein generation</span>
        </div>
        <RouterLink type="button" class="btn btn-primary btn-md mt-2" to="/amino-acid-lab">Amino acid lab</RouterLink>
      </div>

      <!-- Drug Discovery Lab -->
      <div class="col-md-3 justify-content-center">
        <div class="mx-auto protein-container" ref="drugDiscoveryProtein"></div>
        <div class="lab-features">
          <span class="tag">Protein-Ligand Binding</span>
          <span class="upcoming-feature">Metabolic processes affected</span>
        </div>
        <RouterLink type="button" class="btn btn-primary btn-md mt-2" to="/drug-target-lab">Drug discovery lab</RouterLink>
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
  border-radius: 10px; /* Adjust the radius as needed */
}
.lab-features li {
  margin-bottom: 5px; /* Space between list items */
  position: center; /* Needed for positioning the pseudo-element */
  padding-left: 25px; /* Space for the checkmark */
  text-align: left;
}

.tag {
  display: inline-block;
  background-color: #007bff; /* Example background color */
  color: white;
  padding: 5px 10px;
  margin: 5px;
  border-radius: 15px; /* Creates the pill shape */
  font-size: 0.8em;
}

.upcoming-feature {
  display: inline-block; /* Example background color */
  color: white;
  padding: 5px 10px;
  margin: 5px;
  border-radius: 15px; /* Creates the pill shape */
  font-size: 0.8em;
  background-color: #494e4a; /* Different background color for upcoming features */
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