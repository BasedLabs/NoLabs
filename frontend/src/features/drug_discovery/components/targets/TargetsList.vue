<template>
    <q-list bordered separator>
        <q-item-label header>Targets</q-item-label>
        <q-item>
            <q-btn push class="full-width" color="green">
                + add ligands to all targets
            </q-btn>
        </q-item>
        <q-expansion-item popup expand-separator :content-inset-level="1" v-for="target in targets" :key="target.targetId"
            :label="target.targetName"  @show="() => getLigandsForTarget(target)">
           <TargetDetail :target-metadata="target"> </TargetDetail>
            <q-btn push class="full-width" color="green">
                + add ligand
            </q-btn>
            <div v-if="target.loadingLigands">
                <q-spinner color="primary" />
            </div>
            <q-list bordered v-else>
                <q-expansion-item v-for="ligand in target.ligands" :key="ligand.ligandId" clickable :label="ligand.ligandName">
                    <LigandDetail :ligand-metadata="ligand"> </LigandDetail>
                </q-expansion-item>
            </q-list>
        </q-expansion-item>
    </q-list>
    <q-page-sticky position="bottom-right" :offset="[100, 50]">
        <q-btn size="lg" round color="green" icon="add">
            <q-tooltip>Add targets to the experiment</q-tooltip>
        </q-btn>
    </q-page-sticky>
</template>

<script>
import TargetDetail from 'src/features/drug_discovery/components/targets/TargetDetail.vue'
import LigandDetail from 'src/features/drug_discovery/components/ligands/LigandDetail.vue'
import { useDrugDiscoveryStore } from 'src/features/drug_discovery/storage';


export default {
    name: 'TargetsList',
    components: {
        TargetDetail,
        LigandDetail
    },
    props: {
        experimentId: {
            type: String,
            required: true
        },
        originalTargets: {
            type: Array,
            required: true
        }
    },
    data() {
        return {
            targets: this.originalTargets
        }
    },
    methods: {
        async getLigandsForTarget(target) {
            if (target.ligands && target.ligands.length > 0) {
                return;
            }
            const store = useDrugDiscoveryStore();
            target.loadingLigands = true;
            try {
                target.ligands = await store.fetchLigandsForTarget(this.experimentId, target.targetId);
            } catch (error) {
                console.error('Error fetching ligands:', error);
            } finally {
                target.loadingLigands = false;
            }
        },
    }
}

</script>