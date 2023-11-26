describe('Main page', () => {
    const inputPdbId = '#proteinFileInput';
    const inputAminoAcidSequenceId = '#inputSequence';
    const submitAminoAcidSequenceId = '#submitInference';
    const localisationImageClass = '.localisation-image';
    const foldingSelector = '#viewport div';

    it('Amino Acid lab opened', () => {
        cy.on('uncaught:exception', () => false);
        cy.visit('amino-acid-lab');
    });

    it('Inference result is visible', () => {
        cy.visit('amino-acid-lab');
        cy.get('.btn.btn-outline-success.add-experiments-button', {timeout: 90000}).first().click();
        cy.get('.bi.bi-check2.btn.btn.btn-outline-success', {timeout: 90000}).first().click();
        cy.get('#submitInference', {timeout: 90000}).first().click();
        cy.get('#mithochondria-list-item', {timeout: 90000}).should('exist');
    });
})