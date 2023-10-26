describe('Main page', () => {
    const inputPdbId = '#proteinFileInput';
    const inputAminoAcidSequenceId = '#inputSequence';
    const submitAminoAcidSequenceId = '#submitInference';
    const localisationImageClass = '.localisation-image';
    const foldingSelector = '#viewport div';

    it('Amino Acid lab opened', () => {
        cy.visit('http://localhost:5174/amino-acid-lab');
    });

    it('Inference result is visible', () => {
        cy.visit('http://localhost:5174/amino-acid-lab');
        cy.get('.btn.btn-outline-success.add-experiments-button', {timeout: 90000}).click();
        cy.get('.bi.bi-check2.btn.btn.btn-outline-success', {timeout: 90000}).click();
        cy.get('#submitInference', {timeout: 90000}).click();
        cy.intercept('**').as('all');
        cy.wait('@all', {timeout: 90000}).its('response.statusCode').should('be.oneOf', [200]);
        cy.get('#mithochondria-list-item', {timeout: 90000}).should('exist');
    });
})