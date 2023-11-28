<template>
  <div class="text-center" v-if="this.$route.path === '/'">
    <div class="labs-flex-container">
      <!-- Amino Acid Lab -->
      <div class="lab-container">
        <div class="protein-container" ref="aminoAcidProtein"></div>
        <RouterLink type="button" class="btn btn-primary btn-md" to="/amino-acid-lab">Amino acid lab</RouterLink>
      </div>
      <!-- Drug Discovery Lab -->
      <div class="lab-container">
        <div class="protein-container" ref="drugDiscoveryProtein"></div>
        <RouterLink type="button" class="btn btn-primary btn-md" to="/drug-target-lab">Drug discovery lab</RouterLink>
      </div>
    </div>
  </div>
  <RouterView />
</template>

<script>
import * as $3Dmol from '3dmol';

export default {
  watch: {
    '$route' (to, from) {
      // This will be called also when the route changes and the component is reused
      if (to.path === '/') {
        this.initProteins();
      }
    }
  },
  methods: {
    initProteins() {
      this.loadProteinModel(this.$refs.aminoAcidProtein, 'exampleProteins/protein.pdb', true);
      this.loadProteinModel(this.$refs.drugDiscoveryProtein, 'exampleProteins/protein_ligand.pdb', false);
    },
    loadProteinModel(element, pdbPath, isCartoon) {
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
          viewer.setStyle({hetflag: true}, {stick: {radius: 0.6, colorscheme: 'greenCarbon'}});
          }
          viewer.zoomTo();
          this.rotateModel(viewer);
          viewer.render();
        });
    },
    rotateModel(viewer) {
      function rotate() {
        viewer.rotate(0.3, {x:0, y:1, z:0}); // Adjust rotation axis and speed as needed
        viewer.render();
        requestAnimationFrame(rotate);
      }
      rotate();
    }
  }
};
</script>

<style>
.btn {
  /* Your existing styles for buttons */
  border-radius: 10px; /* Adjust the radius as needed */
}

</style>